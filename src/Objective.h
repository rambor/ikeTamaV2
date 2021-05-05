//
// Copyright 2020 Robert P. Rambo
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

#ifndef IKETAMA_OBJECTIVE_H
#define IKETAMA_OBJECTIVE_H
#include <string>
#include <vector>
#include <math.h>
#include <regex>
#include <iostream>
#include "sastools/include/DataBase.h"
#include "sastools/include/PofRData.h"
#include "sastools/include/IofQData.h"
#include "sastools/include/Utils.h"

class Objective {
    // HASH KEY => DATA
    std::map <std::string, DataBase *> datasets;
    std::vector<std::string> keys;
    std::vector<float> weights;

public:
    Objective();

    // copy constructor
    Objective(const Objective &toCopy) {

        keys.resize(toCopy.keys.size());
        weights.resize(toCopy.weights.size());
        std::copy(toCopy.keys.begin(), toCopy.keys.end(), keys.begin());
        std::copy(toCopy.weights.begin(), toCopy.weights.end(), weights.begin());

        for(auto & elements : toCopy.datasets){
            if (elements.second->getIsPr()){
                datasets.emplace(elements.first, elements.second->clone());
            } else {
                datasets.emplace(elements.first, elements.second->clone());
            }
        }
    }

    void swap(Objective & other) noexcept {
        using std::swap;
        std::swap(keys, other.keys);
        std::swap(weights, other.weights);
        std::swap(datasets, other.datasets);
    }

    // copy assignment
    Objective & operator=(const Objective & nodeToCopy) {
        Objective tmp(nodeToCopy);
        tmp.swap(*this);
        return *this;
    }

    ~Objective(){
        datasets.clear();
    }

    void addDataObject(std::string datfile);
    DataBase * getMainDataset();
};


#endif //IKETAMA_OBJECTIVE_H
