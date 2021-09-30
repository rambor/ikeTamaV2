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

#ifndef IKETAMA_LATTICEPOINT_H
#define IKETAMA_LATTICEPOINT_H

#include <math.h>
#include <random>
#include <string>
#include <vector>
#include <sastools/vector3.h>


class LatticePoint {
    int index, bins;

    std::vector<float> amplitudes;
    std::vector<unsigned int> occurences;
    std::vector<float> probabilities;
    float weighted_amplitude;
    int total_amplitudes;

public:
    LatticePoint(int index, int bins);

    void addToCounter(float value);

    float guessAmplitude(float rnd_number) const;

    void resetCounter();

    void updateProbabilities();

    void printProbabilities();

    std::string PDBLineForOutput(vector3 & vec);

    int getIndex(){ return index;}

    float getProbability(int index);

    void resetProbabilities(){
        float prob = 1.0f/(float)total_amplitudes;
        for(auto & pp : probabilities){
            pp = prob;
        }
    }

    float getMostProbableAmplitude() const {
        float temp, max = probabilities[0];
        int selected=0;
        for(int i=1; i<total_amplitudes; i++){
            temp = probabilities[i];
            if (max < temp){
                max = temp;
                selected =i;
            }
        }

        return amplitudes[selected];
    }

    float getMaxProbability() const {
        float temp, max = probabilities[0];
        for(int i=1; i<total_amplitudes; i++){
            temp = probabilities[i];
            if (max < temp){
                max = temp;
            }
        }

        return max;
    }

    void setWeightedAmplitude(){
        weighted_amplitude=0;
        for(int i=0; i<total_amplitudes; i++){
            weighted_amplitude += probabilities[i] * amplitudes[i];
        }
    }

    float getWeightedAmplitude(){
        return weighted_amplitude;
    }

    float getLastAmplitude(){
        return amplitudes[amplitudes.size()-1];
    }

    float CDF(float value);

    std::string getDetails();

};


#endif //IKETAMA_LATTICEPOINT_H
