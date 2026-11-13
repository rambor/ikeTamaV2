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
#include "DensityRefiner.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace densrefine {

static inline double d2(const V3 &a, const V3 &b){ double dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z; return dx*dx+dy*dy+dz*dz; }

unsigned int Data::freeCount() const { return (unsigned int)std::count(isFree.begin(), isFree.end(), 1); }

Data makeData(const std::vector<double> &q, const std::vector<double> &I, const std::vector<double> &sigma, int rebin){
    std::vector<double> qq, II, ss;
    for (size_t i=0; i<q.size() && i<I.size() && i<sigma.size(); i++){
        if (q[i] <= 0 || sigma[i] <= 0 || !std::isfinite(I[i])) continue;
        qq.push_back(q[i]); II.push_back(I[i]); ss.push_back(sigma[i]);
    }
    if (rebin <= 0) rebin = std::max(1, (int)std::lround(qq.size()/150.0));
    Data D;
    for (size_t i = 0; i + rebin <= qq.size(); i += rebin){       // average groups of 'rebin' points
        double sq=0, sI=0, sv=0;
        for (int k=0; k<rebin; k++){ sq += qq[i+k]; sI += II[i+k]; sv += ss[i+k]*ss[i+k]; }
        D.q.push_back(sq/rebin); D.I.push_back(sI/rebin); D.s.push_back(std::sqrt(sv)/rebin);
    }
    for (size_t i=0; i<D.q.size(); i++) D.isFree.push_back((char)(i % 2));
    return D;
}

// HCP lattice (a1 = (a,0,0), a2 = (a/2, a sqrt3/2, 0), c = a sqrt(8/3)) inside a box
static std::vector<V3> hcpLattice(double a, V3 lo, V3 hi){
    std::vector<V3> P;
    const double c = a*std::sqrt(8.0/3.0);
    const V3 a1{a,0,0}, a2{0.5*a, 0.5*std::sqrt(3.0)*a, 0}, b{0.5*a, a/(2*std::sqrt(3.0)), 0.5*c};
    int n1 = (int)std::ceil((hi.x-lo.x)/a) + 4, n2 = (int)std::ceil((hi.y-lo.y)/(0.5*std::sqrt(3.0)*a)) + 4, n3 = (int)std::ceil((hi.z-lo.z)/c) + 4;
    for (int i=-n1; i<=n1; i++) for (int j=-n2; j<=n2; j++) for (int k=-n3; k<=n3; k++){
        for (int s=0; s<2; s++){
            V3 p{i*a1.x + j*a2.x + s*b.x, j*a2.y + s*b.y, k*c + s*b.z};
            if (p.x>=lo.x && p.x<=hi.x && p.y>=lo.y && p.y<=hi.y && p.z>=lo.z && p.z<=hi.z) P.push_back(p);
        }
    }
    return P;
}

// cell hash for neighbour searches (cells of edge 'h')
namespace {
struct CellHash {
    double h; std::unordered_map<long long, std::vector<int>> cells;
    static long long key(long x, long y, long z){ return ((x + 1048576LL) << 42) | ((y + 1048576LL) << 21) | (z + 1048576LL); }
    CellHash(const std::vector<V3> &P, double h_) : h(h_){
        for (int i=0; i<(int)P.size(); i++) cells[key(std::lround(std::floor(P[i].x/h)), std::lround(std::floor(P[i].y/h)), std::lround(std::floor(P[i].z/h)))].push_back(i);
    }
    template<class F> void forNear(const V3 &p, F f) const {
        long cx = std::lround(std::floor(p.x/h)), cy = std::lround(std::floor(p.y/h)), cz = std::lround(std::floor(p.z/h));
        for (long x=cx-1; x<=cx+1; x++) for (long y=cy-1; y<=cy+1; y++) for (long z=cz-1; z<=cz+1; z++){
            auto it = cells.find(key(x,y,z)); if (it == cells.end()) continue;
            for (int j : it->second) f(j);
        }
    }
};
}

Refiner::Refiner(const std::vector<V3> &beadsIn, const Data &data, const Options &options) : opt(options), D(data) {
    if (beadsIn.empty()) throw std::invalid_argument("DensityRefiner: no beads");
    if (opt.massKDa <= 0) throw std::invalid_argument("DensityRefiner: molecular mass (kDa) must be > 0");
    if (D.q.size() < 10) throw std::invalid_argument("DensityRefiner: too few q points");
    if (opt.lambdas.empty()) opt.lambdas = {30.0};
    a = opt.spacing;
    const double Vcell = a*a*a/std::sqrt(2.0);
    const double Vdry = opt.specvol*1.66054*opt.massKDa*1000.0;   // A^3

    // centre the aligned set; outputs are shifted back to the input frame
    centre = {0,0,0};
    for (auto &p : beadsIn){ centre.x+=p.x; centre.y+=p.y; centre.z+=p.z; }
    centre.x/=beadsIn.size(); centre.y/=beadsIn.size(); centre.z/=beadsIn.size();
    std::vector<V3> beads(beadsIn);
    for (auto &p : beads){ p.x-=centre.x; p.y-=centre.y; p.z-=centre.z; }

    // ---- occupancy (KDE) of the aligned set on an HCP lattice
    V3 lo{1e9,1e9,1e9}, hi{-1e9,-1e9,-1e9};
    for (auto &p : beads){ lo.x=std::min(lo.x,p.x); lo.y=std::min(lo.y,p.y); lo.z=std::min(lo.z,p.z); hi.x=std::max(hi.x,p.x); hi.y=std::max(hi.y,p.y); hi.z=std::max(hi.z,p.z); }
    const double pad = 3*a; lo.x-=pad; lo.y-=pad; lo.z-=pad; hi.x+=pad; hi.y+=pad; hi.z+=pad;
    std::vector<V3> grid = hcpLattice(a, lo, hi);
    const double sig = 0.6*a, cut2 = 9*sig*sig;               // KDE kernel, 3 sigma cut-off
    std::vector<double> kde(grid.size(), 0.0);
    {
        CellHash bh(beads, 3*sig);
        #pragma omp parallel for schedule(dynamic, 64)
        for (int g=0; g<(int)grid.size(); g++){
            // same terms as a full loop over the beads, summed in bead order
            std::vector<int> near; bh.forNear(grid[g], [&](int j){ near.push_back(j); });
            std::sort(near.begin(), near.end());
            double s=0; for (int j : near){ double r2=d2(grid[g],beads[j]); if (r2 < cut2) s += std::exp(-r2/(2*sig*sig)); }
            kde[g]=s;
        }
    }
    const double kmax = *std::max_element(kde.begin(), kde.end());
    // support: kde > supportFraction of max, dilated by 'dilate' lattice shells (room for the hydration layer)
    std::vector<char> in(grid.size(), 0);
    for (size_t g=0; g<grid.size(); g++) if (kde[g] > opt.supportFraction*kmax) in[g]=1;
    {
        CellHash gh(grid, 1.01*a);
        for (int pass=0; pass<opt.dilate; pass++){
            std::vector<char> nin = in;
            for (size_t g=0; g<grid.size(); g++) if (in[g]) gh.forNear(grid[g], [&](int h){ if (!in[h] && d2(grid[g],grid[h]) < 1.01*a*a) nin[h]=1; });
            in.swap(nin);
        }
    }
    std::vector<double> occ;
    for (size_t g=0; g<grid.size(); g++) if (in[g]){ pts.push_back(grid[g]); occ.push_back(kde[g]/kmax); }
    const int M = (int)pts.size();
    // prior density: occupancy scaled so that the total equals the dry volume (capped at 1)
    target = Vdry/Vcell;
    double lo_s=0, hi_s=100;
    for (int it=0; it<60; it++){ double s=0.5*(lo_s+hi_s), t=0; for (double o : occ) t += std::min(1.0, s*o); if (t < target) lo_s=s; else hi_s=s; }
    rho0.resize(M); for (int i=0;i<M;i++) rho0[i] = std::min(1.0, lo_s*occ[i]);
    sumRho0 = 0; for (double r : rho0) sumRho0 += r;
    // neighbours (12 at most), listed in index order
    nbr.resize(M); nnbr.assign(M, 0);
    {
        CellHash ph(pts, 1.01*a);
        for (int i=0;i<M;i++){
            std::vector<int> near; ph.forNear(pts[i], [&](int j){ if (j!=i && d2(pts[i],pts[j]) < 1.01*a*a) near.push_back(j); });
            std::sort(near.begin(), near.end());
            for (int j : near) if (nnbr[i] < 12) nbr[i][nnbr[i]++] = j;
        }
    }
    // distance shells: key = round(3 d^2/a^2)
    long maxKey = 0;
    #pragma omp parallel for reduction(max:maxKey) schedule(dynamic, 16)
    for (int i=0;i<M;i++) for (int j=i+1;j<M;j++) maxKey = std::max(maxKey, std::lround(3.0*d2(pts[i],pts[j])/(a*a)));
    std::vector<char> present(maxKey+1, 0); present[0]=1;
    for (int i=0;i<M;i++) for (int j=i+1;j<M;j++) present[std::lround(3.0*d2(pts[i],pts[j])/(a*a))] = 1;
    shellOfKey.assign(maxKey+1, -1);
    for (long k=0;k<=maxKey;k++) if (present[k]){ shellOfKey[k] = (int)shellD.size(); shellD.push_back(a*std::sqrt(k/3.0)); }
    const int K = (int)shellD.size();
    if (K > 65535) throw std::runtime_error("DensityRefiner: too many distance shells");
    if ((size_t)M*M*2 < (size_t)1<<31){            // pair -> shell table (2 bytes per pair) when it fits in 2 GB
        SM.resize((size_t)M*M);
        #pragma omp parallel for schedule(static)
        for (int i=0;i<M;i++) for (int j=0;j<M;j++) SM[(size_t)i*M+j] = (uint16_t)shellOfKey[std::lround(3.0*d2(pts[i],pts[j])/(a*a))];
        Msz = M;
    }
    // shell x q table: bead form factor F(q) of a sphere with the cell volume
    nq = (int)D.q.size();
    const double R = std::cbrt(3*Vcell/(4*M_PI));
    T.assign((size_t)K*nq, 0.0);
    for (int k=0;k<K;k++) for (int iq=0; iq<nq; iq++){
        double q = D.q[iq], x = q*R, F = (x < 1e-6) ? 1.0 : 3*(std::sin(x) - x*std::cos(x))/(x*x*x);
        double qd = q*shellD[k];
        T[(size_t)k*nq + iq] = F*F*((qd < 1e-9) ? 1.0 : std::sin(qd)/qd);
    }
    if (opt.verbose){
        printf("  densrefine: %zu beads | lattice %.3f A, support %d points, %d distance shells | prior total %.0f cells (dry volume %.0f A^3)\n",
               beads.size(), a, M, K, sumRho0, Vdry);
        printf("  densrefine: %d q points (q %.4f to %.4f), %d work / %d free\n", nq, D.q.front(), D.q.back(), nq - (int)D.freeCount(), (int)D.freeCount());
    }
}

void Refiner::setPrTarget(const std::vector<double> &probabilityPerBin, const std::function<int(double)> &binOfDistance){
    nbins = (int)probabilityPerBin.size();
    double norm = 0; for (double v : probabilityPerBin) norm += v;
    if (nbins < 3 || norm <= 0) throw std::invalid_argument("DensityRefiner: bad P(r) target");
    prTarget.clear(); for (double v : probabilityPerBin) prTarget.push_back(v/norm);
    shellBin.assign(shellD.size(), nbins);
    for (size_t k=1; k<shellD.size(); k++){ int b = binOfDistance(shellD[k]); shellBin[k] = std::max(0, std::min(b, nbins)); }   // nbins = beyond dmax
    wP = opt.prWeight;
    if (opt.verbose) printf("  densrefine: P(r) target %d bins, weight %.1f\n", nbins, wP);
}

// IKETAMA PofRData::convertToBin
int binIketama(double distance, double binWidth){
    float ratio = (float)distance/(float)binWidth;
    auto floored = (float)std::floor(ratio);
    if (std::abs(ratio-floored) <= 1000*FLT_EPSILON) return (floored > 0) ? (int)floored - 1 : 0;
    return (int)floored;
}

int Refiner::shellIndex(int i, int j) const {
    if (Msz) return SM[(size_t)i*Msz + j];
    return shellOfKey[std::lround(3.0*d2(pts[i], pts[j])/(a*a))];
}

double Refiner::ampOf(const State &st, int i) const { return st.rho[i] + eta*(1.0 - st.rho[i])*(1.0 - std::exp(-kappa*st.csum[i])); }

void Refiner::computeAmplitudes(State &st) const {
    const size_t M = pts.size();
    st.csum.resize(M); st.amp.resize(M);
    for (size_t i=0; i<M; i++){ double c=0; for (int k=0; k<nnbr[i]; k++) c += st.rho[nbr[i][k]]; st.csum[i]=c; }
    for (size_t i=0; i<M; i++) st.amp[i] = ampOf(st, (int)i);
}

void Refiner::fullShellSums(State &st) const {
    st.S.assign(shellD.size(), 0.0);
    const int M = (int)pts.size();
    for (int i=0; i<M; i++){
        st.S[0] += st.amp[i]*st.amp[i];
        for (int j=i+1; j<M; j++) st.S[shellIndex(i,j)] += 2.0*st.amp[i]*st.amp[j];
    }
}

void Refiner::icalc(const std::vector<double> &Sv, std::vector<double> &out) const {
    out.assign(nq, 0.0);
    for (size_t k=0; k<Sv.size(); k++){
        const double sk = Sv[k]; if (sk == 0.0) continue;
        const double *t = &T[k*nq];
        for (int iq=0; iq<nq; iq++) out[iq] += sk*t[iq];
    }
}

// scale + background on the work set; chi2 work, chi2 free, Durbin-Watson of the work residuals
void Refiner::score(const std::vector<double> &Ic, double &chiW, double &chiF, double &dw, double *cOut, double *bOut) const {
    double sww=0, swx=0, swy=0, swxx=0, swxy=0;
    for (int i=0; i<nq; i++){ if (D.isFree[i]) continue; double w=1.0/(D.s[i]*D.s[i]); sww+=w; swx+=w*Ic[i]; swy+=w*D.I[i]; swxx+=w*Ic[i]*Ic[i]; swxy+=w*Ic[i]*D.I[i]; }
    double det = sww*swxx - swx*swx;
    double c = (sww*swxy - swx*swy)/det, b = (swxx*swy - swx*swxy)/det;
    double sW=0, sF=0; int nW=0, nF=0; double num=0, den=0, prev=0; bool first=true;
    for (int i=0; i<nq; i++){
        double r = (c*Ic[i] + b - D.I[i])/D.s[i];
        if (D.isFree[i]){ sF += r*r; nF++; }
        else { sW += r*r; nW++; den += r*r; if (!first) num += (r-prev)*(r-prev); prev=r; first=false; }
    }
    chiW = sW/std::max(1, nW-2); chiF = sF/std::max(1, nF); dw = den>0 ? num/den : 0;
    if (cOut) *cOut = c;
    if (bOut) *bOut = b;
}

std::vector<double> Refiner::modelPr(const std::vector<double> &Sv) const {
    std::vector<double> P(nbins + 1, 0.0); double tot = 0;
    for (size_t k=1; k<Sv.size(); k++){ P[shellBin[k]] += Sv[k]; tot += Sv[k]; }
    if (tot > 0) for (auto &v : P) v /= tot;
    P.resize(nbins);
    return P;
}

double Refiner::prKL(const std::vector<double> &Sv) const {
    if (prTarget.empty()) return 0.0;
    std::vector<double> P(nbins + 1, 0.0); double tot = 0;
    for (size_t k=1; k<Sv.size(); k++){ P[shellBin[k]] += Sv[k]; tot += Sv[k]; }   // bin nbins = beyond dmax
    if (tot <= 0) return 1e3;
    double kl = 0;
    for (int b=0; b<nbins; b++){ if (prTarget[b] <= 0) continue; double q = P[b]/tot; kl += prTarget[b]*std::log(prTarget[b]/std::max(q, 1e-12)); }
    return kl;
}

double Refiner::rg(const std::vector<double> &rho) const {
    double w=0, x=0, y=0, z=0, r2=0;
    for (size_t i=0; i<pts.size(); i++){ double a_=rho[i]; w+=a_; x+=a_*pts[i].x; y+=a_*pts[i].y; z+=a_*pts[i].z; r2+=a_*(pts[i].x*pts[i].x+pts[i].y*pts[i].y+pts[i].z*pts[i].z); }
    if (w <= 0) return 0;
    x/=w; y/=w; z/=w; return std::sqrt(std::max(0.0, r2/w - x*x - y*y - z*z));
}

std::vector<double> Refiner::modelIntensity(const std::vector<double> &rho) const {
    State st; st.rho = rho; computeAmplitudes(st); fullShellSums(st);
    std::vector<double> Ic; icalc(st.S, Ic);
    double cw, cf, dw, c, b; score(Ic, cw, cf, dw, &c, &b);
    for (auto &v : Ic) v = c*v + b;
    return Ic;
}

double Refiner::selectionThreshold() const {
    double mn = 1e300; for (auto &r : lambdaResults) mn = std::min(mn, r.chiF);
    return mn*(1.0 + 2.0*std::sqrt(2.0/std::max(1u, D.freeCount())));
}

// ------------------------------------------------------------------ snapshots (JSON lines, for animation)
void Refiner::writeSnapshotHeader(FILE *f) const {
    fprintf(f, "{\"type\":\"header\",\"spacing\":%.4f,\"pts\":[", a);
    for (size_t i=0;i<pts.size();i++) fprintf(f, "%s[%.3f,%.3f,%.3f]", i?",":"", pts[i].x+centre.x, pts[i].y+centre.y, pts[i].z+centre.z);   // input frame, as the maps
    fprintf(f, "],\"rho0\":["); for (size_t i=0;i<pts.size();i++) fprintf(f, "%s%.4f", i?",":"", rho0[i]);
    fprintf(f, "],\"q\":["); for (int i=0;i<nq;i++) fprintf(f, "%s%.6e", i?",":"", D.q[i]);
    fprintf(f, "],\"I\":["); for (int i=0;i<nq;i++) fprintf(f, "%s%.6e", i?",":"", D.I[i]);
    fprintf(f, "],\"s\":["); for (int i=0;i<nq;i++) fprintf(f, "%s%.6e", i?",":"", D.s[i]);
    fprintf(f, "],\"free\":["); for (int i=0;i<nq;i++) fprintf(f, "%s%d", i?",":"", (int)D.isFree[i]);
    fprintf(f, "],\"eta\":%.3f,\"nbins\":%d,\"prTarget\":[", eta, nbins);
    for (int b=0;b<nbins;b++) fprintf(f, "%s%.6e", b?",":"", prTarget[b]);
    fprintf(f, "]}\n");
}

void Refiner::writeSnapshot(FILE *f, const State &st, long step, long total, double temperature, double lambda) const {
    std::vector<double> I2; icalc(st.S, I2); double cw, cf, dw, c, b; score(I2, cw, cf, dw, &c, &b);
    fprintf(f, "{\"type\":\"state\",\"step\":%ld,\"total\":%ld,\"T\":%.4e,\"lambda\":%.3f,\"chiW\":%.5f,\"chiF\":%.5f,\"dw\":%.4f,\"kl\":%.4e,\"rg\":%.3f,\"rho\":[",
            step, total, temperature, lambda, cw, cf, dw, prKL(st.S), rg(st.rho));
    for (size_t i=0;i<st.rho.size();i++) fprintf(f, "%s%.3f", i?",":"", st.rho[i]);
    fprintf(f, "],\"Icalc\":["); for (int i=0;i<nq;i++) fprintf(f, "%s%.6e", i?",":"", c*I2[i]+b);
    if (!prTarget.empty()){
        std::vector<double> P = modelPr(st.S);
        fprintf(f, "],\"pr\":["); for (int k=0;k<nbins;k++) fprintf(f, "%s%.6e", k?",":"", P[k]);
    }
    fprintf(f, "]}\n");
}

// ------------------------------------------------------------------ one annealing run at fixed lambda
LambdaResult Refiner::anneal(double lam) const {
    const int M = (int)pts.size(), K = (int)shellD.size();
    const double mu = opt.mu;
    State st; st.rho = rho0; computeAmplitudes(st); fullShellSums(st);
    std::vector<double> Ic, Icand(nq), dS(K), oldAmp; std::vector<int> changed; changed.reserve(13); oldAmp.reserve(13);
    icalc(st.S, Ic);
    double cw, cf, dw; score(Ic, cw, cf, dw);
    std::mt19937 gen(opt.seed); std::uniform_int_distribution<int> pick(0, M-1); std::normal_distribution<double> gauss(0,1); std::uniform_real_distribution<double> uni(0,1);
    double priorSq=0, sumRho=0; for (int i=0;i<M;i++){ double d=st.rho[i]-rho0[i]; priorSq+=d*d; sumRho+=st.rho[i]; }
    auto energyWith = [&](double chiW, double pSq, double sR, double kl){ double v=sR/sumRho0-1.0; return chiW + lam*(pSq/M)/0.01 + mu*v*v + wP*kl; };
    double klCur = prKL(st.S);
    double E = energyWith(cw, priorSq, sumRho, klCur);
    // starting temperature: median |dE| of a few random moves
    const long steps = opt.sweeps*(long)M;
    double Tt = 0;
    {
        std::vector<double> dEs;
        for (int t=0; t<=20; t++){
            int i=pick(gen); double nr=std::min(1.0,std::max(0.0,st.rho[i]+0.3*gauss(gen)));
            double keep=st.rho[i]; st.rho[i]=nr; computeAmplitudes(st); fullShellSums(st); icalc(st.S, Icand); double a1,b1,c1; score(Icand,a1,b1,c1);
            double pS=priorSq-(keep-rho0[i])*(keep-rho0[i])+(nr-rho0[i])*(nr-rho0[i]);
            dEs.push_back(std::fabs(energyWith(a1,pS,sumRho-keep+nr,prKL(st.S))-E)); st.rho[i]=keep;
        }
        computeAmplitudes(st); fullShellSums(st);
        std::sort(dEs.begin(), dEs.end()); Tt = std::max(1e-6, dEs[dEs.size()/2]);
    }
    const double Tstart = Tt, Tend = Tt*1e-3;

    FILE *snap = nullptr;
    if (opt.snapshots > 0 && !opt.snapshotPrefix.empty()){
        char name[64]; snprintf(name, sizeof(name), "_lambda%g_snapshots.jsonl", lam);
        snap = fopen((opt.snapshotPrefix + name).c_str(), "w");
        if (snap) writeSnapshotHeader(snap);
    }
    const long snapEvery = snap ? std::max(1L, steps/opt.snapshots) : 0;
    if (snap) writeSnapshot(snap, st, 0, steps, Tstart, lam);

    // the state kept is the one with the lowest chi2 free (+ P(r) term when a P(r) target is set)
    LambdaResult best; best.lambda = lam; best.chiW = cw; best.chiF = cf; best.dw = dw; best.kl = klCur; best.rho = st.rho;
    double bestSel = cf + wP*klCur;
    long accepted = 0; double step = 0.3;
    const long report = std::max(1L, steps/10);
    for (long s=0; s<steps; s++){
        if (snapEvery && s > 0 && s % snapEvery == 0){ computeAmplitudes(st); fullShellSums(st); writeSnapshot(snap, st, s, steps, Tstart*std::pow(Tend/Tstart, (double)s/steps), lam); }
        const double Tc = Tstart*std::pow(Tend/Tstart, (double)s/steps);
        const int i = pick(gen);
        const double oldR = st.rho[i];
        double nr = oldR + step*gauss(gen); nr = std::min(1.0, std::max(0.0, nr));
        if (nr == oldR) continue;
        // rho_i changes -> amplitude of i, and the neighbour sums (and amplitudes) of its neighbours
        changed.clear(); oldAmp.clear(); std::fill(dS.begin(), dS.end(), 0.0);
        st.rho[i] = nr;
        auto applyAmp = [&](int p, double newA){
            double old = st.amp[p], dlt = newA - old; if (std::fabs(dlt) < 1e-14) return;
            changed.push_back(p); oldAmp.push_back(old);
            dS[0] += newA*newA - old*old;
            if (Msz){ const uint16_t *row = &SM[(size_t)p*Msz]; for (int j=0;j<M;j++){ if (j==p) continue; dS[row[j]] += 2.0*dlt*st.amp[j]; } }
            else for (int j=0;j<M;j++){ if (j==p) continue; dS[shellIndex(p,j)] += 2.0*dlt*st.amp[j]; }
            st.amp[p] = newA;
        };
        for (int k=0;k<nnbr[i];k++) st.csum[nbr[i][k]] += nr - oldR;
        applyAmp(i, ampOf(st, i));
        for (int k=0;k<nnbr[i];k++){ int j=nbr[i][k]; applyAmp(j, ampOf(st, j)); }
        for (int k=0;k<K;k++) dS[k] += st.S[k];
        icalc(dS, Icand);
        double nW, nF, nDW; score(Icand, nW, nF, nDW);
        const double nPS = priorSq - (oldR-rho0[i])*(oldR-rho0[i]) + (nr-rho0[i])*(nr-rho0[i]);
        const double nSR = sumRho - oldR + nr;
        const double nKL = prKL(dS);
        const double nE = energyWith(nW, nPS, nSR, nKL);
        if (nE <= E || uni(gen) < std::exp(-(nE-E)/Tc)){
            st.S.swap(dS); E=nE; priorSq=nPS; sumRho=nSR; accepted++;
            if (nF + wP*nKL < bestSel){ bestSel=nF+wP*nKL; best.chiF=nF; best.kl=nKL; best.rho=st.rho; best.chiW=nW; best.dw=nDW; }
        } else {
            st.rho[i] = oldR;
            for (int k=0;k<nnbr[i];k++) st.csum[nbr[i][k]] -= nr - oldR;
            for (int k=(int)changed.size()-1; k>=0; k--) st.amp[changed[k]] = oldAmp[k];
        }
        if ((s+1) % report == 0){
            double rate = (double)accepted/report; accepted = 0;
            step = std::min(0.5, std::max(0.02, step*(rate > 0.3 ? 1.2 : 0.8)));
            computeAmplitudes(st); fullShellSums(st);                 // resynchronise (floating-point drift)
            if (opt.verbose){
                icalc(st.S, Ic); double a1,b1,c1; score(Ic,a1,b1,c1);
                #pragma omp critical(densrefine_print)
                { printf("    lambda %6.1f %3ld%%  T %.2e  chi2 work %.3f free %.3f  DW %.2f  KL %.2e  Rg %.1f  (best free %.3f)\n",
                         lam, (s+1)*100/steps, Tc, a1, b1, c1, prKL(st.S), rg(st.rho), best.chiF); fflush(stdout); }
            }
        }
    }
    if (snap){ State b; b.rho = best.rho; computeAmplitudes(b); fullShellSums(b); writeSnapshot(snap, b, steps, steps, Tend, lam); fclose(snap); }
    best.rg = rg(best.rho);
    return best;
}

void Refiner::run(){
    // hydration contrast: scan with the prior density
    State st; st.rho = rho0;
    std::vector<double> Ic;
    double bestEta=0, bestChi=1e300;
    for (int e10=0; e10<=12; e10++){
        eta = 0.05*e10; computeAmplitudes(st); fullShellSums(st); icalc(st.S, Ic);
        double cw, cf, dw; score(Ic, cw, cf, dw);
        if (cw < bestChi){ bestChi=cw; bestEta=eta; }
    }
    eta = bestEta;
    computeAmplitudes(st); fullShellSums(st); icalc(st.S, Ic);
    score(Ic, chiW0, chiF0, dw0); kl0 = prKL(st.S);
    if (opt.verbose) printf("  densrefine start (aligned-set occupancy): eta %.2f | chi2 work %.3f free %.3f | DW %.2f | P(r) KL %.2e | Rg %.1f\n",
                            eta, chiW0, chiF0, dw0, kl0, rg(rho0));

    const int L = (int)opt.lambdas.size();
    lambdaResults.assign(L, LambdaResult());
    #pragma omp parallel for schedule(dynamic, 1)
    for (int l=0; l<L; l++) lambdaResults[l] = anneal(opt.lambdas[l]);

    // choice: largest lambda whose chi2 free is within noise of the best
    const double thr = selectionThreshold();
    chosenIndex = 0; double chosenL = -1;
    for (int l=0; l<L; l++) if (lambdaResults[l].chiF <= thr && lambdaResults[l].lambda > chosenL){ chosenL = lambdaResults[l].lambda; chosenIndex = l; }
    if (opt.verbose){
        printf("  densrefine lambda ladder (chi2 free limit %.3f = min x (1 + 2 sqrt(2/%u))):\n", thr, D.freeCount());
        for (int l=0; l<L; l++){ auto &r = lambdaResults[l];
            printf("    lambda %6.1f: chi2 work %.3f free %.3f DW %.2f P(r) KL %.2e Rg %.1f%s\n", r.lambda, r.chiW, r.chiF, r.dw, r.kl, r.rg, (size_t)l == chosenIndex ? "   <- chosen" : ""); }
    }
}

// ------------------------------------------------------------------ output
void Refiner::writeXPLOR(const std::string &fn, const std::vector<double> &rho) const {
    const double h = 0.5*a, sigma = 0.5*a;
    std::vector<V3> P(pts); for (auto &p : P){ p.x+=centre.x; p.y+=centre.y; p.z+=centre.z; }
    double lo[3]={1e9,1e9,1e9}, hi[3]={-1e9,-1e9,-1e9};
    for (auto &p : P){ double c[3]={p.x,p.y,p.z}; for (int k=0;k<3;k++){ lo[k]=std::min(lo[k],c[k]); hi[k]=std::max(hi[k],c[k]); } }
    int mn[3], mx[3];
    for (int k=0;k<3;k++){ mn[k]=(int)std::floor((lo[k]-3*sigma)/h); mx[k]=(int)std::ceil((hi[k]+3*sigma)/h); }
    int nx=mx[0]-mn[0]+1, ny=mx[1]-mn[1]+1, nz=mx[2]-mn[2]+1;
    std::vector<double> G((size_t)nx*ny*nz, 0.0);
    int rc = (int)std::ceil(3*sigma/h);
    for (size_t i=0; i<P.size(); i++){
        if (rho[i] <= 0) continue;
        int ix=(int)std::lround(P[i].x/h)-mn[0], iy=(int)std::lround(P[i].y/h)-mn[1], iz=(int)std::lround(P[i].z/h)-mn[2];
        for (int z=std::max(0,iz-rc); z<=std::min(nz-1,iz+rc); z++) for (int y=std::max(0,iy-rc); y<=std::min(ny-1,iy+rc); y++) for (int x=std::max(0,ix-rc); x<=std::min(nx-1,ix+rc); x++){
            double dx=(x+mn[0])*h-P[i].x, dy=(y+mn[1])*h-P[i].y, dz=(z+mn[2])*h-P[i].z;
            G[((size_t)z*ny + y)*nx + x] += rho[i]*std::exp(-(dx*dx+dy*dy+dz*dz)/(2*sigma*sigma));
        }
    }
    FILE *f = fopen(fn.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + fn);
    fprintf(f, "\n       2 !NTITLE\nREMARK  densrefine: refined density against I(q)\nREMARK  rho on HCP lattice, Gaussian-splatted\n");
    fprintf(f, "%8d%8d%8d%8d%8d%8d%8d%8d%8d\n", nx, mn[0], mx[0], ny, mn[1], mx[1], nz, mn[2], mx[2]);
    fprintf(f, "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\nZYX\n", nx*h, ny*h, nz*h, 90.0, 90.0, 90.0);
    double mean=0, sd=0; for (double v : G) mean += v; mean /= G.size(); for (double v : G) sd += (v-mean)*(v-mean); sd = std::sqrt(sd/G.size());
    for (int z=0; z<nz; z++){
        fprintf(f, "%8d\n", z);
        size_t n=0;
        for (int y=0; y<ny; y++) for (int x=0; x<nx; x++){ fprintf(f, "%12.5E", G[((size_t)z*ny+y)*nx+x]); if (++n % 6 == 0) fprintf(f, "\n"); }
        if (n % 6) fprintf(f, "\n");
    }
    fprintf(f, "   -9999\n%12.4E%12.4E\n", mean, sd);
    fclose(f);
}

void Refiner::writePDB(const std::string &fn, const std::vector<double> &rho, const std::vector<int> &which, const std::string &remark) const {
    FILE *f = fopen(fn.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + fn);
    fprintf(f, "REMARK 265 %s\nREMARK 265                   EDGE RADIUS : %.3f\n", remark.c_str(), 0.5*a);
    int n=1;
    for (int i : which){
        fprintf(f, "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n", n%100000, n%10000,
                pts[i].x+centre.x, pts[i].y+centre.y, pts[i].z+centre.z, 100.0*rho[i]);
        n++;
    }
    fprintf(f, "END\n"); fclose(f);
}

void Refiner::write(const std::string &prefix) const {
    if (lambdaResults.empty()) throw std::logic_error("DensityRefiner::write before run");
    const LambdaResult &best = chosen();
    const int M = (int)pts.size();
    State st; st.rho = best.rho; computeAmplitudes(st); fullShellSums(st);
    std::vector<double> Ic; icalc(st.S, Ic);
    double cw, cf, dw, scale, bg; score(Ic, cw, cf, dw, &scale, &bg);
    {
        FILE *f = fopen((prefix + "_fit.dat").c_str(), "w");
        if (!f) throw std::runtime_error("cannot write " + prefix + "_fit.dat");
        fprintf(f, "# q I_obs sigma I_calc free  (chi2 work %.4f free %.4f DW %.3f eta %.2f lambda %.2f)\n", cw, cf, dw, eta, best.lambda);
        for (int i=0;i<nq;i++) fprintf(f, "%.6e %.6e %.6e %.6e %d\n", D.q[i], D.I[i], D.s[i], scale*Ic[i]+bg, (int)D.isFree[i]);
        fclose(f);
    }
    std::vector<int> pos; for (int i=0;i<M;i++) if (best.rho[i] > 0.01) pos.push_back(i);
    writePDB(prefix + "_density.pdb", best.rho, pos, "refined density: B-factor = 100 x rho");
    std::vector<int> ord(M); for (int i=0;i<M;i++) ord[i]=i;
    std::stable_sort(ord.begin(), ord.end(), [&](int x, int y){ return best.rho[x] > best.rho[y]; });
    double acc=0; std::vector<int> top; for (int i : ord){ if (acc >= target) break; top.push_back(i); acc += 1.0; }
    writePDB(prefix + "_beads.pdb", best.rho, top, "refined bead model: highest-density points, dry volume");
    writeXPLOR(prefix + "_map.xplor", best.rho);
    writeXPLOR(prefix + "_start_map.xplor", rho0);
}

} // namespace densrefine
