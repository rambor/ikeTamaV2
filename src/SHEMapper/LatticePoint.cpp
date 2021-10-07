//
// Created by xos81802 on 13/06/2021.
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


#include <iostream>
#include "LatticePoint.h"


LatticePoint::LatticePoint(int index, int bins) : index(index), bins(bins) {

//    float delta = 1.0f/(float)bins;
//    for (int i=0; i<bins; i++){
//        amplitudes.push_back(i*delta);
//    }

    amplitudes.resize(bins);

    amplitudes[0] = -0.269f;
    amplitudes[1] = 0.6345f;
    amplitudes[2]=1.0f;


    // electron densities - nucleic, protein, lipid, and lipid subtracted water and divided by largest value
//    amplitudes[0] = -0.269f;
//    amplitudes[1] = 0.103f;
//    amplitudes[2] = 0.682543f;
//    amplitudes[3] = 1.0f;

    total_amplitudes = amplitudes.size();
    probabilities.resize(total_amplitudes);
    occurences.resize(total_amplitudes);

    this->resetProbabilities();
    this->setWeightedAmplitude();
}

/*
 * should be a random number between 0 and 1
 */
float LatticePoint::guessAmplitude(float random_number) const {
    float sum = probabilities[0];
    for(int i=1; i<total_amplitudes; i++){ // CDF - cumalative distribution function
        if (random_number <= sum){
            return amplitudes[i-1];
        }
        sum+=probabilities[i];
    }

    return amplitudes[total_amplitudes-1];
    //return CDF(random_number);
}


float LatticePoint::CDF(float prob){
    float sum = probabilities[0];
    for(int i=1; i<total_amplitudes; i++){ // CDF - cumalative distribution function
        if (prob <= sum){
            return amplitudes[i-1];
        }
        sum+=probabilities[i];
    }

    return amplitudes[total_amplitudes-1];
}


void LatticePoint::resetCounter(){

    for(auto & cnt : occurences){
        cnt = 0;
    }
}

void LatticePoint::addToCounter(float value) {
    /*
     * amp
     * 0
     * 0.25
     * 0.5
     * 0.75
     * 1.0
     *
     *     amplitudes[0] = -0.1f;
     *     amplitudes[1] = 0.103f;
     *     amplitudes[2] = 0.682543f;
     *     amplitudes[3] = 1.0f;
     *
     */
    for(int i=0; i<total_amplitudes; i++){
        if (abs(value - amplitudes[i]) < 0.01){
            occurences[i]+=1;
            return;
        }
    }
}


void LatticePoint::updateProbabilities(){

    int sum = 0;
    for(auto & point : occurences){
        sum += point;
    }

    auto invTotal = 1.0f/(float)sum;

    const float updateAlpha=0.67f;
    float oldprob, newprob;

    int index = 0;
    float temp;
    for(auto & point : occurences){

        oldprob = (1.0f - updateAlpha)*probabilities[index];

        if (point > 0){
            newprob = updateAlpha*((float)point*invTotal) + oldprob;
        } else {
            newprob = oldprob;
        }

        probabilities[index] = newprob;
        index++;
    }

}



void LatticePoint::printProbabilities(){

    float best = probabilities[0];
    float sum = 0.0f;

    int count = 0, keep=0;
    for(auto & pp : probabilities){
        std::cout << count << " " << amplitudes[count] << " " << pp << " " << occurences[count] << std::endl;
        sum += pp;
        if (pp > best){
            best = pp;
            keep = count;
        }
        count++;
    }
    std::cout << "LATTICE :: "<< index << " sum " << " " << sum << std::endl;
}


std::string LatticePoint::PDBLineForOutput(vector3 & vec){

    float best = probabilities[0];
    int count = 0, keep=0;
    for(auto & pp : probabilities){
        if (pp > best){
            best = pp;
            keep = count;
        }
        count++;
    }

    printf("%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00%5.2f\n", "ATOM", 1," CA ", "ALA", "A", index, vec.x, vec.y, vec.z, probabilities[keep] );
    return "string";
}

float LatticePoint::getProbability(int index){
    return probabilities[index];
}


std::string LatticePoint::getDetails(){

    std::string tempHeader = "";
    char buffer[80];

    tempHeader.append("Probability and Amplitude Details");
    for(std::size_t i=0; i < amplitudes.size(); i++){
        std::snprintf(buffer, 80, "%5i %.5E %.5E \n", ((int)i+1), amplitudes[i], probabilities[i]);
        tempHeader.append(buffer);
    }

    return tempHeader;
}