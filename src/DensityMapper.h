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

// #include <boost/math/special_functions/spherical_harmonic.hpp

class DensityMapper {

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

    void setDebyeFactors( std::vector<Datum> & workingSet);

    void populateICalcOpenMP(unsigned int total_q, unsigned int totalHCPInUse, std::vector<float> & iCalc, float * squared_amplitudes, float * db);
    void populateICalc(unsigned int total, std::vector<float> & pICalc, float * squared_amplitudes);

    float populateDensities(std::vector<float> &amplitudesT, std::vector<float> &densities);
    void populateDensitiesOMP(unsigned int totalHCPInUse, std::vector<Neighbors> neighbors, std::vector<float> & amplitudesT, std::vector<float> & densities_at_HCP);

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

};


#endif //IKETAMA_DENSITYMAPPER_H
