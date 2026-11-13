//
// Created by xos81802 on 10/06/2021.
//
// Copyright 2021 Robert P. Rambo
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

#ifndef IKETAMA_DENSITYMAPPER_H
#define IKETAMA_DENSITYMAPPER_H

#include <omp.h>
#include <math.h>
#include <random>
#include <sastools/Datum.h>
#include <Eigen/Core>
#include <sastools/simde-no-tests-master/simde-arch.h>
#include <sastools/simde-no-tests-master/simde-common.h>

#include <string>
#include <vector>
#include <sastools/vector3.h>
#include "SHEMapper/Shell.h"
#include "boost/math/special_functions/bessel.hpp"
#include "SHEMapper/LatticePoint.h"
#include "PointSetModel.h"
#include "SHEMapper/Neighbors.h"
#include <iostream>
#include <sastools/PDBModel.h>
#include <sastools/utils.h>
#include <cstdlib>
#include <cstdint>
#include <cfloat>
#include "DensityRefiner.h"

class PofRData;

// #include <boost/math/special_functions/spherical_harmonic.hpp

class DensityMapper {
    friend struct DensityMapperTestAccess;   // unit tests

    struct Trial{
        double value;
        float dw;
        float counter;
        float scale;
        std::vector<float> model_amplitudes;
        std::vector<float> residuals;
        Trial() = default;

        Trial(double val, float dw, float cntr, float scale, std::vector<float>  & vec, std::vector<float>  & res) : value(val), dw(dw), counter(cntr), scale(scale) {
            //indices = std::vector<unsigned int>(vec);
            model_amplitudes = std::move(vec);
            residuals = std::move(res);
        }

        bool operator<(const Trial & a) const{
            return value < a.value;
        }
    };

    struct Density{
        unsigned int index;
        float density;
        Density() = default;
        Density(unsigned int index, float val) : index(index), density(val) {
        }
    };

    float sampling_frequency;
    float delta_r;
    float cutoff;
    float kmax;
    float fwhm_sigma;
    float ylm_size;
    float within_limit, upper_limit;

    std::vector<float> rvalues;
    std::vector<double> cos_thetas;
    std::vector<double> phis;
    std::vector<float> bessels;
    std::vector<double> y_lm_real;
    std::vector<double> y_lm_imag;

    std::vector<Neighbors> neighborhoods;
    //std::vector<float> debye_factors;
    float * debye_factors;
    //aligned_vector debye_factors;

    PointSetModel modelHCPSamplingDensityPoints;//, inputBaseModel;

    std::vector<Shell> shells;
    std::vector<LatticePoint> lattice_points;

    std::vector<vector3> centered_coordinates;

    std::vector<float> amplitudes;

    std::map<unsigned int, std::set<unsigned int >> neighboringHCPPointsOfModelLattice;
    std::vector<unsigned int> keptHCPSamplingPoints;

    std::vector<bool> useIt;

    int total_shells, lmax, bessel_size, qvalues_size, total_centered_coordinates;

    float dmin_supremum, dmin_infimum, stdev_dmin, average_dmin;
    float max_neighbors;

    void createSphericalFibonacciLattice();

    /*
     * I(q) from shell sums. The sampling points lie on an HCP lattice, so their pair distances take a small number of
     * values (shells). I(q) = sum_k S_k sinc(q r_k) with S_k = sum of rho_i*rho_j over the pairs in shell k (2x for i<j,
     * shell 0 = the diagonal). Replaces the (K(K+1)/2 x Nq) Debye-factor table: O(K^2/2 + shells*Nq) per evaluation and
     * a few MB of memory instead of a table of K^2/2*Nq floats (copied once per thread).
     */
    std::vector<uint32_t> pairShell;      // shell of pair (i<j), upper triangle row by row
    std::vector<float> shellRadius;       // r_k
    std::vector<float> shellSinc;         // sinc(q r_k), laid out [q][k]
    unsigned int totalPairShells = 0, shellTableQ = 0;
    unsigned long rngSeed = 0;
    std::vector<int> keptIndexOfHCPPoint;  // HCP universe index -> position in keptHCPSamplingPoints (-1 if not kept)
    std::string inputFilename;             // aligned set (EDGE RADIUS remark gives the bead lattice for refineModelDensity)

public:
    /*
     * given a set of aligned structures,
     */
    DensityMapper(std::string filename, float qmax, float sampling_frequency);

    // create radial, spherical lattice for sampling real space densities
    // r, theta, phi
    // need delta_r, theta and phi will come from

    /*
     *
     * each round is guess a set of amplitudes
     * score them and update probability model
     *
     */

    void setModel();

    int getTotalShells(){ return total_shells;}

    void setBessels(const std::vector<float> &qvalues);

    void calculateDensityCoefficientAtLMR();

    float calculateAmplitudeAtQ(int q_index, int l_index, int m_index);
    float calculateIntensityAtQ(int q_index);

    void setAmplitudes(std::vector<float> &amplitudes);

    void refineModel(int max_rounds, float topPercent, int models_per_round,
                     std::vector<Datum> & workingSet,
                     std::vector<Datum> & workingSetSmoothed);

    void printLatticePointsInfo();

    int getTotalCenteredCoordinates(){ return total_centered_coordinates;}

    std::string headerParametersForXPLOR(int &na, int &nb, int &nc, float grid_spacing);
    std::string headerParametersForXPLORFlipped(int & na, int & nb, int & nc, float grid_spacing);

    void createXPLORMap(std::string name);

    void tester(std::string pdbfile,  std::vector<Datum> & workingSet);

    void trimWhiteSpace(std::string &text);

    void createHCPGrid();

    void setDebyeFactors( std::vector<Datum> & workingSet);   // legacy (table of K^2/2 x Nq), not used by the refinements

    // shell tables for the q values of workingSet (call after createHCPGrid)
    void setShellTables(std::vector<Datum> & workingSet);
    // I(q) of the densities at the kept sampling points; thread-safe, shellSums is caller-owned scratch
    void calculateICalcFromDensities(const float * densities, float * iCalc, std::vector<double> & shellSums) const;
    unsigned int getTotalPairShells() const { return totalPairShells; }

    // reproducible runs: 0 (default) seeds from std::random_device as before
    void setRandomSeed(unsigned long seed){ rngSeed = seed; }
    unsigned int seedFor(unsigned int stream){
        if (rngSeed != 0) return (unsigned int)(rngSeed + 7919u*stream);
        std::random_device rd; return rd();
    }

    void populateICalcOpenMP(unsigned int total_q, unsigned int totalHCPInUse, std::vector<float> & iCalc, float * squared_amplitudes, float * db);
    void populateICalc(unsigned int total, std::vector<float> & pICalc, float * squared_amplitudes);

    float populateDensities(std::vector<float> &amplitudesT, std::vector<float> &densities);
    void populateDensitiesOMP(unsigned int totalHCPInUse, std::vector<Neighbors> & neighbors, std::vector<float> & amplitudesT, std::vector<float> & densities_at_HCP);

    float calculateScaleFactor(unsigned int total, float *const pICalc, float *const pSigma, Datum *const pWorkingSet);

    float getChiSquare(unsigned int total, float scale, float *const pICalc, float *const pSigma,
                       Datum *const pWorkingSet, Datum *const pWorkingSetSmoothed, float * const res);

    void writeICalc(unsigned int total, float scale,
                    float *const pICalc,
                    Datum *const pWorkingSet,
                    Datum * const pWorkingSetSmoothed,
                    std::string name);

    float calculateDurbinWatson(unsigned int total, const float * const residuals);

    void openMP();

    void writeLatticePoints(std::string name);

    float calculate_rmsd(std::vector<float> & priors);

    unsigned int probability_count(float limit);
    unsigned int getTotalKeptHCPSamplingPoints(){return keptHCPSamplingPoints.size();}

    void refineModelASA(unsigned int highTempRounds, std::vector<Datum> &workingSet,
                   std::vector<Datum> &workingSetSmoothed);

    void updateASATemp(unsigned int index, unsigned int evalMax, float acceptRate, double &temp, double &inv_temp);

    void updateDensities(std::vector<unsigned int> &indicesInUse, std::vector<float> &amplitudesT,
                         std::vector<float> &densities_at_HCP, float *squared_amplitudes);

    const std::vector<LatticePoint> &getLatticePoints() const {
        return lattice_points;
    }


    void updateICalc(std::vector<Density> & priorHCPdensities,
                     const unsigned int total_data,
                     std::vector<float> & i_calc,
                     float * prior_hcp_electron_densities,
                     float * hcp_electron_densities);

    void getIndicesOfHCPSamplingGrid(std::vector<unsigned int> & indices_selected,
                                     std::vector<unsigned int> & hcpSamplingPointsToUpdate);

    void populateSquaredAmplitudes(unsigned int totalHCPInUse, std::vector<float> &densities_at_HCP,
                                   float *squared_amplitudes);

    void printParameters(std::vector<float> &temp,
                         std::vector<double> &accept,
                         std::vector<float> &score,
                         std::vector<float> & chis, std::vector<int> & flips);

    double f1(double x, double y);

    void refineModelOPENMP(int max_threads, int max_rounds, float topPercent, int models_per_round,
                           std::vector<Datum> &workingSet, std::vector<Datum> &workingSetSmoothed);

    float convolutionFunction(float length);

    void setMaxNeighbors();

    void setLatticePointNeighbors();

    bool checkConnectivity(std::set<unsigned int> &selectedIndices, std::vector<LatticePoint> & omp_lattice_points);


    std::vector<LatticePoint> & getLatticePoints(){ return lattice_points;}

    /*
     * Density refinement against I(q) (and P(r)) starting from the occupancy of the aligned set; replaces the
     * cross-entropy search (refineModel / refineModelOPENMP) with Metropolis annealing of a continuous density on the
     * bead lattice (see DensityRefiner.h). A ladder of prior weights (options.lambdas) is run in parallel and the
     * largest weight whose chi2(free) is within noise of the best is kept.
     *
     *   workingSet   I(q) to fit (q, I, sigma); rebinned to ~150 points, alternate points held out as the free set
     *   pofrData     optional (nullptr): P(r) restraint, KL(P(r)_data || P(r)_model) on the Shannon bins
     *   options      massKDa is required; spacing <= 0 takes 2 x EDGE RADIUS of the aligned set
     *   prefix       writes <prefix>_map.xplor, _start_map.xplor, _density.pdb, _beads.pdb, _fit.dat
     *
     * Coordinates are in the centred frame of the aligned set (as CenteredInputModel).
     */
    densrefine::LambdaResult refineModelDensity(std::vector<Datum> & workingSet, PofRData * pofrData,
                                                densrefine::Options options, const std::string & prefix = "refined_density");

    static float readEdgeRadius(const std::string & pdbfile);

};


#endif //IKETAMA_DENSITYMAPPER_H
