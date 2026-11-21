//
// Copyright 2026 Robert P. Rambo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
// DensityRefiner: refine a density on an HCP lattice against I(q) (and optionally P(r)), starting from the
// occupancy (KDE) of an aligned set of bead models. Used by DensityMapper::refineModelDensity and by the
// stand-alone densrefine program.
//
//   density     rho_i in [0,1] on HCP lattice points (spacing a); support = KDE of the aligned set above a
//               threshold, dilated by a few lattice shells
//   amplitude   a_i = rho_i + eta (1 - rho_i)(1 - exp(-1.5 * sum_{nbr j} rho_j))   (hydration layer, contrast eta)
//   I(q)        = c * F(q)^2 * sum_k S_k sinc(q d_k) + b,  S_k = sum of a_i a_j over the pairs in distance shell k.
//               HCP distances are discrete (3 d^2/a^2 is an integer), so a single-point change updates the S_k in
//               O(points) and I(q) in O(shells x q); no Debye-factor table. Lattice shells closer than
//               Options::shellMerge (0.2 A) are merged at their pair-weighted mean distance (error in sinc ~ (q dr)^2/24,
//               < 1e-4 at q = 0.25), which cuts the shell count about 3x.
//   data        I(q) is averaged down to Options::pointsPerShannon points per Shannon channel (dmax qmax / pi) before
//               fitting; the cost of every move is proportional to (shells x q points).
//   energy      chi2(work) + lambda * mean((rho - rho_KDE)^2)/0.01 + mu * (sum rho / sum rho_KDE - 1)^2
//               + w * KL(P(r)_data || P(r)_model)          (last term only when a P(r) is given)
//               c and b by weighted least squares at every evaluation; alternate q points are held out (free set).
//   search      Metropolis single-point moves, temperature falls 1000x geometrically. For each lambda the state with
//               the lowest chi2(free) (+ w*KL) is kept.
//   lambda      a ladder of values is run (in parallel); the chosen lambda is the LARGEST one whose chi2(free) is within
//               statistical noise of the best:  chi2_free <= min chi2_free * (1 + 2 sqrt(2 / N_free)).
//               Small lambda lets the density drift into diffuse shapes that fit I(q) equally well (seen for dsDNA);
//               this picks the model closest to the aligned set that the data cannot distinguish from the best fit.
//
#ifndef IKETAMA_DENSITYREFINER_H
#define IKETAMA_DENSITYREFINER_H

#include <array>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace densrefine {

struct V3 { double x, y, z; };

struct Data {                         // I(q) after rebinning; isFree marks the held-out points
    std::vector<double> q, I, s;
    std::vector<char> isFree;
    unsigned int freeCount() const;
};
// average groups of 'rebin' points (0 = choose so that about 150 points remain); alternate points are held out
// rebin factor that leaves about perShannon points per Shannon channel ((qmax - qmin) dmax / pi, at least 3 channels)
int rebinForShannon(const std::vector<double> &q, double dmax, int perShannon);
Data makeData(const std::vector<double> &q, const std::vector<double> &I, const std::vector<double> &sigma, int rebin = 0);
// bin of a distance, IKETAMA PofRData::convertToBin convention (lower < d <= upper)
int binIketama(double distance, double binWidth);

struct Options {
    double spacing = 5.023;           // HCP lattice spacing (A), normally the bead lattice of the models
    double massKDa = 0;               // molecular mass (kDa); with specvol sets the total density (dry volume)
    double specvol = 0.73;            // partial specific volume (cm^3/g): 0.73 protein, ~0.55 DNA/RNA
    double mu = 10.0;                 // total-density restraint
    double prWeight = 1000.0;         // weight of the P(r) term (used only when a P(r) target is set)
    double supportFraction = 0.02;    // support: KDE > fraction of its maximum ...
    double shellMerge = 0.2;          // merge distance shells closer than this (A); 0 = exact lattice shells
    int pointsPerShannon = 6;         // I(q) points per Shannon channel after rebinning (half are held out); 0 = ~150 points
    int dilate = 2;                   // ... dilated by this many lattice shells
    long sweeps = 100;                // moves per support point
    unsigned int seed = 1;
    int snapshots = 0;                // > 0: write this many states per lambda to <prefix>_lambda<L>_snapshots.jsonl
    std::string snapshotPrefix;
    std::vector<double> lambdas{0, 1, 3, 10, 30, 100, 300};   // one value = no selection
    bool verbose = true;
};

struct LambdaResult {
    double lambda = 0, chiW = 0, chiF = 0, dw = 0, kl = 0, rg = 0;
    std::vector<double> rho;          // kept state
};

class Refiner {
    friend struct DensityRefinerTestAccess;   // unit tests
public:
    Refiner(const std::vector<V3> &beads, const Data &data, const Options &options);

    // P(r) target: probability per bin (bins 0..n-1) and the bin of a distance (>= n means beyond dmax).
    // With PofRData: prob[b] = getProbabilityPerBin(b), binOf = convertToBin.
    void setPrTarget(const std::vector<double> &probabilityPerBin, const std::function<int(double)> &binOfDistance);

    void run();                                        // all lambdas, then the choice
    void write(const std::string &prefix) const;       // _map.xplor _start_map.xplor _density.pdb _beads.pdb _fit.dat

    // results
    const std::vector<LambdaResult> & results() const { return lambdaResults; }
    const LambdaResult & chosen() const { return lambdaResults[chosenIndex]; }
    double startChiW() const { return chiW0; }
    double startChiF() const { return chiF0; }
    double startDW() const { return dw0; }
    double startKL() const { return kl0; }
    double hydrationContrast() const { return eta; }
    unsigned int supportPoints() const { return (unsigned int)pts.size(); }
    unsigned int distanceShells() const { return (unsigned int)shellD.size(); }
    const std::vector<V3> & points() const { return pts; }
    const std::vector<double> & priorDensity() const { return rho0; }
    // model I(q) (scaled, with background) of a density on the support points
    std::vector<double> modelIntensity(const std::vector<double> &rho) const;
    double selectionThreshold() const;                 // chi2_free limit used for the choice

private:
    struct State {                                     // everything one annealing run changes
        std::vector<double> rho, amp, csum, S;
    };
    Options opt;
    Data D;
    V3 centre{0,0,0};                                  // input frame = lattice frame + centre
    double a = 0, eta = 0, kappa = 1.5, sumRho0 = 0, target = 0;
    std::vector<V3> pts;
    std::vector<std::array<int,12>> nbr; std::vector<int> nnbr;
    std::vector<double> rho0;
    std::vector<int> shellOfKey; std::vector<uint16_t> SM; int Msz = 0;
    std::vector<double> shellD; int latticeShells = 0;
    std::vector<double> T;                             // T[k*nq + iq] = F(q)^2 sinc(q d_k)
    int nq = 0;
    std::vector<double> prTarget; std::vector<int> shellBin; int nbins = 0; double wP = 0;
    double chiW0 = 0, chiF0 = 0, dw0 = 0, kl0 = 0;
    std::vector<LambdaResult> lambdaResults; size_t chosenIndex = 0;

    int shellIndex(int i, int j) const;
    double ampOf(const State &st, int i) const;
    void computeAmplitudes(State &st) const;
    void fullShellSums(State &st) const;
    void icalc(const std::vector<double> &Sv, std::vector<double> &out) const;
    void score(const std::vector<double> &Ic, double &chiW, double &chiF, double &dw, double *c = nullptr, double *b = nullptr) const;
    double prKL(const std::vector<double> &Sv) const;
    std::vector<double> modelPr(const std::vector<double> &Sv) const;
    double rg(const std::vector<double> &rho) const;
    LambdaResult anneal(double lambda) const;
    void writeSnapshotHeader(FILE *f) const;
    void writeSnapshot(FILE *f, const State &st, long step, long total, double temperature, double lambda) const;
    void writeXPLOR(const std::string &fn, const std::vector<double> &rho) const;
    void writePDB(const std::string &fn, const std::vector<double> &rho, const std::vector<int> &which, const std::string &remark) const;
};

} // namespace densrefine

#endif //IKETAMA_DENSITYREFINER_H
