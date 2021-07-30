//
// Created by xos81802 on 18/06/2021.
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

#ifndef IKETAMA_SHELL_H
#define IKETAMA_SHELL_H


#include <sastools/include/vector3.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <set>

//extern "C" {
//   void plmbar_(double *    plm,
//             int *       lmax,
//             double *    z,
//             int *       csphase,
//             int *       cnorm
//   );
//};


class Shell {

    int total_coordinates_in_shell;
    std::vector<vector3> points; // cartesian coordinates
    std::vector<bool> useIt;
    std::vector<float> thetas;
    std::vector<float> phis;
    std::vector<float> distances;
    std::vector<double> plms;

    float radius, area_per_point, dmin;
public:
    float getRadius() const;

private:

    int ylm_size, lmax;
    std::vector<float> y_lm_real;
    std::vector<float> y_lm_imag;
    std::vector<float> p_lm_real;
    std::vector<float> p_lm_imag;

    std::map<unsigned int long, std::vector<unsigned int long >> neighbors;
    void createShell();

public:

    Shell(float radius_shell, float area_per_point);

    float getEpsilon();

    void setModel(std::vector<vector3> & centered_coordinates, float cutoff, float sigma);

    int getPointsPerShell(){ return total_coordinates_in_shell;}

    void populateSphericalHarmonics(int lmax);

    void updateDensityCoefficients(std::vector<vector3> & centered_coordinates, std::vector<float> & amplitudes, float sigma);

    float getPLMReal(int l_index, int m_index);
    float getPLMImag(int l_index, int m_index);
};


#endif //IKETAMA_SHELL_H
