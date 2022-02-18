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
#include <set>


class LatticePoint {
    int index, bins;

    std::vector<float> amplitudes;
    std::vector<unsigned int> occurences;
    std::vector<float> probabilities;
    std::set<unsigned int> neighbors;

    float weighted_amplitude;
    unsigned int total_amplitudes;
public:
    unsigned int getTotalAmplitudes() const;

private:
    unsigned int assigned_index;

public:
    LatticePoint(int index, int bins);

    // copy constructor
    LatticePoint(const LatticePoint & dat) : index(dat.index), bins(dat.bins) {
        this->amplitudes = dat.amplitudes;
        this->occurences = dat.occurences;
        this->probabilities = dat.probabilities;
        this->neighbors = dat.neighbors;
        this->weighted_amplitude = dat.weighted_amplitude;
        this->total_amplitudes = dat.total_amplitudes;
        this->assigned_index = dat.assigned_index;
    }

    LatticePoint & operator=(const LatticePoint & dataToCopy) {
        if (this == &dataToCopy)
            return *this;

        LatticePoint tmp(dataToCopy); // make a copy
        tmp.swap(*this);
        return *this;
    }

    /**
     * Rule of 3.5, define copy, destructor and assignment operator
     * @param other
     */
    void swap(LatticePoint & other) {
        other.index = index;
        other.bins = bins;

        other.amplitudes = std::move(amplitudes);
        other.occurences = std::move(occurences);
        other.probabilities = std::move(probabilities);
        other.neighbors = std::move(neighbors);

        other.weighted_amplitude = weighted_amplitude;
        other.total_amplitudes = total_amplitudes;
        other.assigned_index = assigned_index;
    }

    void addToCounter(float value);

    float guessAmplitude(float rnd_number) const;

    float getAmplitudeByIndex(unsigned int index) const { return amplitudes[index];}

    void resetCounter();

    void updateProbabilities(float updateAlpha=0.67f);

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
            if (temp > max){
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

    float getWeightedAmplitude() const {
        return weighted_amplitude;
    }

    float getLastAmplitude(){
        return amplitudes[amplitudes.size()-1];
    }

    float CDF(float value);

    std::string getDetails();

    void setRandomIndexForConvergence(){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<unsigned int> randomIndex(0,total_amplitudes-1);
        assigned_index = randomIndex(gen);
    }

    float getAssignedAmplitude() const {
        return probabilities[assigned_index];
    }

    void updateOccurence(unsigned int indexToUpdate){
        occurences[indexToUpdate] += 1;
    }

    void calculateProbabilitiesFromOccurrences();

    void addNeighbor(unsigned int index);

    unsigned int getTotalNeighbors(){ return neighbors.size();}

    std::set<unsigned int> & getNeighbors(){ return neighbors;}

    bool isNeighbor(unsigned int index){
        return (neighbors.find(index) != neighbors.end());
    }

};


#endif //IKETAMA_LATTICEPOINT_H
