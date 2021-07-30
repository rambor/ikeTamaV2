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

#include <math.h>
#include <random>
#include <sastools/include/Datum.h>
#include <Eigen/Core>
#include <sastools/include/simde-no-tests-master/simde-arch.h>
#include <sastools/include/simde-no-tests-master/simde-common.h>
#include <sastools/include/simde-no-tests-master/x86/avx.h>
//#include <immintrin.h>

#ifdef __APPLE__ //Mac OSX has a different name for the header file
    #include <OpenCL/opencl.h>
#else
    #include <CL/cl.h>
#endif

#include <string>
#include <vector>
#include <sastools/include/vector3.h>
#include "SHEMapper/Shell.h"
#include "boost/math/special_functions/bessel.hpp"
#include "SHEMapper/LatticePoint.h"
#include "PointSetModel.h"
#include "SHEMapper/Neighbors.h"
#include <iostream>
#include <sastools/include/PDBModel.h>
#include <sastools/include/utils.h>
#include <cstdlib>

// #include <boost/math/special_functions/spherical_harmonic.hpp

class DensityMapper {

    struct Trial{
        double value;
        float counter;
        float scale;
        std::vector<float> model_amplitudes;
        std::vector<float> residuals;
        Trial() = default;

        Trial(double val, float cntr, float scale, std::vector<float>  & vec, std::vector<float>  & res) : value(val), counter(cntr), scale(scale) {
            //indices = std::vector<unsigned int>(vec);
            model_amplitudes = std::move(vec);
            residuals = std::move(res);
        }

        bool operator<(const Trial & a) const{
            return value < a.value;
        }
    };

    float sampling_frequency;
    float delta_r;
    float cutoff;
    float kmax;
    float fwhm_sigma;
    float ylm_size;

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

    PointSetModel modelDensityHCP;

    std::vector<Shell> shells;
    std::vector<LatticePoint> lattice_points;

    std::vector<vector3> centered_coordinates;
    std::vector<vector3> grid_for_map;
    std::vector<float> amplitudes;

    std::map<unsigned int long, std::vector<unsigned int long >> neighbors;
    std::vector<unsigned int> keptHCPLatticePoints;

    std::vector<bool> useIt;

    int total_shells, lmax, bessel_size, qvalues_size, total_centered_coordinates;

    float dmin_supremum, dmin_infimum, stdev_dmin, average_dmin;

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

    void refineModel(int max_rounds, float topPercent, int models_per_round,  std::vector<Datum> & workingSet);

    void createDensityMapGrid();

    void printLatticePointsInfo();

    int getTotalCenteredCoordinates(){ return total_centered_coordinates;}

    std::string headerParametersForXPLOR(int &na, int &nb, int &nc);
    std::string headerParametersForXPLORFlipped(int & na, int & nb, int & nc);

    void createXPLORMap(std::string name);

    void tester(std::string pdbfile,  std::vector<Datum> & workingSet);

    void trimWhiteSpace(std::string &text);

    void createHCPGrid();

    void setDebyeFactors( std::vector<Datum> & workingSet);

    void populateICalc(unsigned int total, std::vector<float> & pICalc, float * squared_amplitudes);

    float populateDensities(std::vector<float> &amplitudesT, std::vector<float> &densities, float * squared_amplitudes);

    float calculateScaleFactor(unsigned int total, float *const pICalc, float *const pSigma, Datum *const pWorkingSet);

    float getChiSquare(unsigned int total, float scale, float *const pICalc, float *const pSigma,
                       Datum *const pWorkingSet, float * const res);

    void writeICalc(unsigned int total, float scale, float *const pICalc, Datum *const pWorkingSet, std::string name);

    float calculateDurbinWatson(unsigned int total, float * const residuals);

    int getTotal_centered_coordinates() const;

//    void openMP();

};


#endif //IKETAMA_DENSITYMAPPER_H
