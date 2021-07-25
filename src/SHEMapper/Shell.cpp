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


#include <set>
#include "Shell.h"

Shell::Shell(float radius, float area_per_point): radius(radius), area_per_point(area_per_point) {


    total_coordinates_in_shell = (radius > 0) ? (int)std::ceil((radius*radius*4.0*M_PI)/area_per_point) : 1;

    //std::cout << " SHELL " << radius << " " << total_coordinates_in_shell << std::endl;
    createShell();
}


void Shell::createShell() {

    points.resize(total_coordinates_in_shell);
    thetas.resize(total_coordinates_in_shell);
    phis.resize(total_coordinates_in_shell);

    vector3 * shell_vec;

    auto inv_golden_ratio = (float)((2.0f*M_PI)/((1.0f + sqrtf(5.0))*0.5f));

    float epsilon = getEpsilon(), phi_calc;

    for(int j=0; j < total_coordinates_in_shell; j++){

        shell_vec = &points[j];
        double theta = j*inv_golden_ratio;
        double phi = acos(1.0d - 2.0d*(j+epsilon)/(total_coordinates_in_shell - 1 + 2.0d*epsilon));

        double sin_phi = sin(phi);
        shell_vec->x = radius*(float)(cos(theta)*sin_phi);
        shell_vec->y = radius*(float)(sin(theta)*sin_phi);
        shell_vec->z = radius*(float)(cos(phi));

        // reassign angles to match coordinate system
        phi_calc = (shell_vec->x == 0 || shell_vec->x == 0.0d) ? 0.0d : atan2f(shell_vec->y,shell_vec->x);
        // theta = phi; //acosf((shell_vec->z/radius)); // acosf => range 0 to PI
        thetas[j] = (float)phi;
        phis[j] = phi_calc;
    }


    dmin = 0;
    float least, dmin_count = 0.0f;

    for(int i=0; i< (total_coordinates_in_shell-1); i++){
        shell_vec = &points[i];
        least = 1000;
        for(int j=(i+1); j<total_coordinates_in_shell; j++) {

            float temp = (*shell_vec- points[j]).length();

            if (temp < least){
                least = temp;
            }
        }

        dmin += least;
        dmin_count += 1.0f;
    }

    dmin = (dmin > 0) ? dmin/dmin_count : 0;
    //std::cout << radius <<  " shell dmin " << dmin << std::endl;

}


float Shell::getEpsilon(){

    if (total_coordinates_in_shell >= 600000){
        return 214;
    } else if (total_coordinates_in_shell >= 400000){
        return 75;
    } else if (total_coordinates_in_shell >= 11000){
        return 27;
    } else if (total_coordinates_in_shell >= 890){
        return 10;
    } else if (total_coordinates_in_shell >= 177){
        return 3.33;
    } else if (total_coordinates_in_shell >= 24){
        return 1.33;
    }

    return 0.33;
}

/**
 *
 * @param centered_coordinates pointset model of supporting structure
 * @param dmin_lattice sets size of neighborhood to consider in calculationn
 * @param sigma
 */
void Shell::setModel(std::vector<vector3> & centered_coordinates, float dmin_lattice, float sigma) {

    float dis, cutoff =  (dmin_lattice+dmin)*0.5f;
    float inv_sigma = 1.0d/(2.0*sigma*sigma);

    useIt.resize(total_coordinates_in_shell);

    // which shell coordinates are not within cut-off centered_coordinates
    unsigned int long total_in_coords = centered_coordinates.size();
    int count=0;

    /*
     * find all the points in the model that are in the neighborhood of the shell lattice points
     * for each shell point determine closest model points and hole inn a container
     */
    for(unsigned int long i=0; i<total_coordinates_in_shell; i++){

        vector3 * shell_vec = &points[i];

        bool flagit = true;
        for(unsigned int long j=0; j<total_in_coords; j++){
            vector3 * coor_vec = &centered_coordinates[j];
            dis = (*(shell_vec) - *(coor_vec)).length();

            if (dis < cutoff){ // sets up which points to include as neighbors in the calculation
                useIt[i] = true;
                count++;
                flagit = false;
                break;
            }
        }

        if (flagit){
            useIt[i] = false;
        }
    }

    // for each point in lattice, determine its neighbors in the centered coordinates
    for(unsigned int long i=0; i<total_coordinates_in_shell; i++){
        if (useIt[i]){
            vector3 * shell_vec = &points[i];
            for(unsigned int long j=0; j<total_in_coords; j++){
                vector3 * coor_vec = &centered_coordinates[j];
                dis = (*(shell_vec) - *(coor_vec)).length();

                if ( dis < cutoff){ // cutoff should be more specific to
                    auto it = neighbors.find(i);
                    if (it == neighbors.end()){
                        auto temp = neighbors.emplace(i, std::vector<unsigned int long>() );
                        (*temp.first).second.push_back(j);
                    } else {
                        it->second.push_back(j);
                    }

                    distances.push_back(exp(-dis*dis*inv_sigma));
                }
            }
        }
    }

    int sum=0;
    for (auto & pts : neighbors){
        auto & neighborhood = pts.second;
        sum += neighborhood.size();
//        std::cout << radius << " " << pts.first << " " << neighborhood.size() << " SHELL " << total_coordinates_in_shell << std::endl;
    }
    //std::cout <<  "SHELL " << radius << " IN USE " << neighbors.size() << " <NGHOOD SIZE> " << ((float)sum/(float)neighbors.size()) << std::endl;
}

void Shell::populateSphericalHarmonics(int lmax){

    this->lmax = lmax;
    int csphase = -1, cnorm = 1;
    int plmSize = ((lmax+1)*(lmax+2))/2 + 1; // (lmax+1)*(lmax+2)/2, add one to componsate fo

    ylm_size = (lmax+1)*(lmax+1)*total_coordinates_in_shell;
    y_lm_real.resize((unsigned long)ylm_size);
    y_lm_imag.resize((unsigned long)ylm_size);
    p_lm_real.resize((unsigned long)(lmax+1)*(lmax+1));
    p_lm_imag.resize((unsigned long)(lmax+1)*(lmax+1));

    plms.resize((unsigned long)(total_coordinates_in_shell*(plmSize-1)));
    std::vector<double> plm((unsigned long)plmSize);
    double * const ptrPLM = plms.data();

    unsigned int plmCount=0;

    double costheta;

//    costheta = cos(M_PI/180.0*45);
//    plmbar_(&plm[0], &lmax, &costheta, &csphase, &cnorm);
//    std::cout << " LMAX " << lmax << std::endl;
//    for (int j=1; j < plmSize; j++){
//        std::cout << j << " " << (float)plm[j] << std::endl;
//    }

    for(unsigned int long i=0; i<total_coordinates_in_shell; i++){

        if (useIt[i]){
            costheta = cos(thetas[i]);
            // lmax
            plmbar_(&plm[0], &lmax, &costheta, &csphase, &cnorm); //FORTRAN SHTOOLS SHE
            for (int j=1; j < plmSize; j++){
                ptrPLM[plmCount] = (float)plm[j];
                plmCount++;
            }

        } else {
            for (int j=1; j < plmSize; j++){
                ptrPLM[plmCount] = 0.0d;
                plmCount++;
            }
        }
    }

    unsigned int index_lm = 0;
    double factor = 0.5d * (sqrt(1.0/M_PI));
    float sinp, cosp;
    double result;
    float phiValue;
    float * const ptrR = y_lm_real.data();
    float * const ptrI = y_lm_imag.data();

    for(int l=0; l <= lmax; l++) { //  y_lm's
        for(int m=-l; m <= l; m++) { //assemble sphherical harmonics at constant l,m for all atoms

            int lm_index = (l*(l+1))/2 + abs(m);

            for (int i=0 ; i < total_coordinates_in_shell; i++) { // can't assume numAtoms is even
                // includes negative angles for m < 0
                phiValue = (-m*phis[i]); // Y*_lm is (-1)^m*Y_(-m, l)
                result = ptrPLM[(i*(plmSize-1) + lm_index)] * factor;
                __sincosf(phiValue, &sinp, &cosp);
                ptrR[ index_lm ] = (float) (cosp * result);
                ptrI[ index_lm ] = (float) (sinp * result);
                index_lm++;
            }
        }
    }
}




void Shell::updateDensityCoefficients(std::vector<vector3> & centered_coordinates, std::vector<float> & amplitudes, float sigma){

    float density;

    float inv_sqrt = 1.0d/(sigma*sqrt(2.0*M_PI));
    float inv_sigma = 1.0d/(2.0*sigma*sigma);

    int index;
    float * const ptrR = y_lm_real.data();
    float * const ptrI = y_lm_imag.data();

    float * const pPlmR = p_lm_real.data();
    float * const pPlmI = p_lm_imag.data();

    float * const pAmps = amplitudes.data();
    float * const pExp = distances.data();

    float real_sum, imag_sum;
    int p_lm_count = 0;
    int dis_index;

    for(int l=0; l <= lmax; l++){

        int start_at = l*l*total_coordinates_in_shell;

        for(int m=-l; m<=l; m++){

            int m_spot = start_at + (l + m)*total_coordinates_in_shell;

            real_sum = 0.0;
            imag_sum = 0.0;
            dis_index=0;

            for(auto & it : neighbors){
                density = 0.0;
                //vector3 & shell_vec  = points[it.first];
                for (auto & jit : it.second) {  // go over neighborhood (intersection of spherical lattice with model)
                    // could be sped up due to distance calculations being fixed
//                    vector3 *coor_vec = &centered_coordinates[(unsigned int) jit];
//                    float dis = (shell_vec - *(coor_vec)).length();
                    //density += amplitudes[(unsigned int) jit] * inv_sqrt * exp(-dis * dis * inv_sigma);
                    //density += amplitudes[(unsigned int) jit] * inv_sqrt * exp(-dis * dis * inv_sigma);
                    density += pAmps[(unsigned int) jit] * inv_sqrt * pExp[dis_index];
                    //std::cout << pExp[dis_index] << " == " << exp(-dis * dis * inv_sigma) << std::endl;
                    dis_index++;
                }
                index = m_spot+it.first;

                real_sum += ptrR[ index ]*density;
                imag_sum += ptrI[ index ]*density;
            }

//            for(unsigned int long i=0; i<total_coordinates_in_shell; i++){
//                if (useIt[i]) {
//                    auto it = neighbors.find(i);
//                    density = 0.0;
//                    for (auto & jit : it->second) {  // go over neighborhood (intersection of spherical lattice with model)
//                        // could be sped up due to distance calculations being fixed
//                        //vector3 *coor_vec = &centered_coordinates[(unsigned int) jit];
//                        //dis = (*(shell_vec) - *(coor_vec)).length();
//                        //density += amplitudes[(unsigned int) jit] * inv_sqrt * exp(-dis * dis * inv_sigma);
//                        //density += amplitudes[(unsigned int) jit] * inv_sqrt * exp(-dis * dis * inv_sigma);
//                        density += pAmps[(unsigned int) jit] * inv_sqrt * pdis[dis_index];
//                        //std::cout << distances[dis_index] << " == " << exp(-dis * dis * inv_sigma) << std::endl;
//                        dis_index++;
//                    }
//
//                    // given coordinate index i, l and get harmonic
//                    index = m_spot+i;
//
//                    real_sum += ptrR[ index ]*density;
//                    imag_sum += ptrI[ index ]*density;
//                }
//            }
            pPlmR[p_lm_count] = real_sum;
            pPlmI[p_lm_count] = imag_sum;
            p_lm_count++;
        }
    }

    //std::cout << "max density " << density << std::endl;
}

float Shell::getRadius() const {
    return radius;
}

float Shell::getPLMReal(int l_index, int m_index) {
    return p_lm_real[l_index*l_index + (l_index+m_index)];
}

float Shell::getPLMImag(int l_index, int m_index) {
    return p_lm_imag[l_index*l_index + (l_index+m_index)];
}
