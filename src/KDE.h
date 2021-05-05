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

#ifndef IKETAMA_KDE_H
#define IKETAMA_KDE_H

#include "../thirdparty/Eigen/Core"
#include "../thirdparty/Eigen/Geometry"
#include "../thirdparty/power_sasa.h"

#include <string>
#include <vector>
#include <algorithm>
#include <sastools/include/vector3.h>
#include <sastools/include/PDBModel.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <iostream>
#include <cstdio>
#include <utility>
#include <limits>
#include <fstream>
#include "PointSetModel.h"
#include "Anneal.h"

class PointSetModel;

struct ProbabilityBead{ double prob; unsigned int index;
    bool operator<(const ProbabilityBead & a) const{
        return prob < a.prob;
    }
};


class KDE {

    struct Trial{ double score; std::vector<unsigned int> indices;};

    bool priorSet = false;
    std::string prior_name;
    unsigned int minN=INT16_MAX, maxN=0, kdeSUM;
    std::string listofPDBFiles;
    std::vector<vector3> coordinates;
    std::vector<vector3> prior;
    std::vector<vector3> centered_prior;
    std::vector<vector3> centered_coordinates;
    std::vector<std::string> pdbfiles;
    std::map<unsigned int, std::vector<unsigned int> > bead_to_vec_map;
    std::map<std::string, float > kernel_mask_mapping;
    std::map<std::string, float > kernel_mapping;
    vector3 centering_vec;
    unsigned int totalCoordinates, totalFiles, trialSize, topNSize;
    float bead_radius, topNpercent;
    float minx =  std::numeric_limits<float>::max();
    float miny =  std::numeric_limits<float>::max();
    float minz =  std::numeric_limits<float>::max();
    float maxx =  -std::numeric_limits<float>::max();
    float maxy =  -std::numeric_limits<float>::max();
    float maxz =  -std::numeric_limits<float>::max();
    float cminx, cminy, cminz, cmaxx, cmaxy, cmaxz;
    float kdeCount, grid_spacing;
    float bandwidth, delta_r = 1.8;
    float average_dmin, dmin_supremum, dmin_infimum, stdev_dmin;
    float solventContrast = 0.334;
    float proteinContrast = 0.414;
    float nucleicacidContrast = 0.621;

    double remappedAverage, remappedStDev;

    float remapToNewHCPLattice(PointSetModel * pModel, float cutoff);
    float remapToNewHCPLatticeSym(PointSetModel * pModel, float cutoff);
    void createMaskOnCenteredObject(std::vector<unsigned int> & mask, float limit);
    std::string headerParametersForXPLOR(int & na, int & nb, int & nc, float & adjustedCMINX, float & adjustedCMINY, float & adjustedCMINZ);
    std::string headerParametersForXPLORSYM(int & na, int & nb, int & nc, float & adjustedCMINX, float & adjustedCMINY, float & adjustedCMINZ, std::string sym);
    std::string mirrorTempHeader;
    void logger(std::string description, std::string value);
    std::string formatNumber(float number, int decimals);

    unsigned int estimateRounds(unsigned int lowerBound, unsigned int totalInSet);
public:

    KDE(std::string filename, float bead_radius, float topNpercent);

    std::string getFileExt(const std::string &s);

    bool checkKDEFile(std::string basic_string);
    bool createKDE(PointSetModel * pModel);

    void extractCoordinates();

    float getCenteredMaxZ() {return std::abs(cminz) > std::abs(cmaxz) ? std::abs(cminz) : std::abs(cmaxz);}
    float getCenteredMaxX() {return std::abs(cminx) > std::abs(cmaxx) ? std::abs(cminx) : std::abs(cmaxx);}
    float getCenteredMaxY() {return std::abs(cminy) > std::abs(cmaxy) ? std::abs(cminy) : std::abs(cmaxy);}

    float kernel(float value);

    float map_refine(PointSetModel *pModel, PofRData *pData, std::string outname);
    float map_refineSym(PointSetModel *pModel, PofRData *pData, std::string outname);

    void generateSeedFromInputFile(PointSetModel *pModel, PofRData *pData, std::string outname);

    void writeCoordinatesToFile(std::string name, std::vector<vector3> & coords) const;

    float getMinx(){return minx;}
    float getMiny(){return miny;}
    float getMinz(){return minz;}
    float getMaxx(){return maxx;}
    float getMaxy(){return maxy;}
    float getMaxz(){return maxz;}

    void add_prior(std::string);
    float getAverageDmin(){ return average_dmin;}
    float getDminStDev(){ return stdev_dmin; }
    float getSupremumDmin(){ return dmin_supremum;}
    float surfaceVolumeCalculation(unsigned int workingLimit, std::vector<unsigned int> &mask, PointSetModel *pModel);


    unsigned int numberOfContactsFromSet(std::set<unsigned int> *beads_in_use,
                                                 PointSetModel *pModel,
                                                 unsigned int const selectedIndex);
};


#endif //IKETAMA_KDE_H
