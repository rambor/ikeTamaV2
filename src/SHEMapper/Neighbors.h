//
// Created by xos81802 on 05/07/2021.
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

#ifndef IKETAMA_NEIGHBORS_H
#define IKETAMA_NEIGHBORS_H

#include <vector>

class Neighbors {

    unsigned int hcp_lattice_index, kept_index;

    std::vector<unsigned int> model_indices;
    std::vector<float> distances;

public:
    Neighbors(unsigned int hcp_lattice_index);

    // copy constructor
    Neighbors(const Neighbors &n2){
        this->hcp_lattice_index = n2.hcp_lattice_index;
        this->kept_index = n2.kept_index;

        for(auto & it : n2.model_indices){
            this->model_indices.push_back(it);
        }

        for(auto & it : n2.distances){
            this->distances.push_back(it);
        }
    }

    void swap(Neighbors & other) noexcept {
        using std::swap;
        other.hcp_lattice_index = hcp_lattice_index;
        other.kept_index = kept_index;
        std::swap(model_indices, other.model_indices);
        std::swap(distances, other.distances);
    }

    Neighbors & operator=(const Neighbors & nodeToCopy) {
        Neighbors tmp(nodeToCopy);
        tmp.swap(*this);
        return *this;
    }

    Neighbors(Neighbors && other) noexcept  {
        hcp_lattice_index = other.hcp_lattice_index;
        kept_index = other.kept_index;
        model_indices = std::move(other.model_indices);
        distances = std::move(other.distances);
    }

    ~Neighbors() = default;

    void add_neighbor(unsigned int index, float distance){
        model_indices.push_back(index);
        distances.push_back(distance);
    }

    void clear_neighbors(){
        model_indices.clear();
        distances.clear();
    }

    const unsigned int * getModelIndices(){ return model_indices.data();}
    const float * getModelDistances(){ return distances.data();}

    const unsigned int getTotalNeighbors() { return model_indices.size();}

    unsigned int getHCPLatticeIndex() { return hcp_lattice_index;}

    void setKeptIndex(unsigned int val){ kept_index = val;}
    unsigned int getKeptIndex(){ return kept_index;}

};


#endif //IKETAMA_NEIGHBORS_H
