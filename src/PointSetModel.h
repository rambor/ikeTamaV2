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

#ifndef IKETAMA_MODEL_H
#define IKETAMA_MODEL_H
#include <string>
#include <vector>
#include <math.h>
#include <regex>
#include <iostream>
#include <sastools/Bead.h>
#include <sastools/PofRData.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <set>
#include "SubUnit.h"
#include "Anneal.h"

class SubUnit;
class Anneal;

class PointSetModel {

    float radius_of_universe, bead_radius, cutOffNeighbor;
    float xaxis, yaxis, zaxis;
    float inv_bead_radius, helicalVolume;
    float beadAverage, beadStDev, cvx_volume;
    float averageNumberOfContactsInModel;

    unsigned int sizeOfNeighborhood = 12;
    unsigned int number_of_beads, total_in_reduced_seed, total_in_working_universe, neighborLimit, total_in_seed; // there will be no bead with this index
    unsigned long int totalDistances, startingWorkingLimit, workingLimit, deadLimit, baseWorkingLimit;
    std::vector<unsigned short int> bins; // linear array (n*(n-1)/2) should this be unsigned short int, bins will never be greater than 100
    float sasa_volume_start=0.0f, sasa_surface_area=0.0f, surface_to_volume=0.0f;
    float radial_limit, bead_volume;

    bool useDirectMethod = false, useCylindricalSearchSpace=false;

    std::string symmetry="C1", helical_point_group_symmetry = "C1";
    float risePerSubUnit, rotationPerSubunit;

    unsigned int numberOfSubUnits=1, numberOfHelicalSubUnits, numberOfHelicalSubUnitsInSingleFilament;
    unsigned int symmetryIndex=1;
    std::string symmetryGroup="";
    std::vector<unsigned int> neighboringSubunits; // subunits that correspond to neighbor of the asymmetry unit (zero order)
    std::vector<Bead> beads;
    std::vector<unsigned int> bead_indices;
    std::vector<unsigned int> neighbors; // 12 * bead_indices size
    std::vector<unsigned int> startingSetVector;
    std::vector<unsigned int> seed_indices;  // contains true model of PDB model converted to lattice model
    std::set<unsigned int> reduced_seed;

    std::vector<float> subUnitAnglesCos;
    std::vector<float> subUnitAnglesSin;

    std::vector<SubUnit> subUnits;

public:

    PointSetModel() = default;
    PointSetModel(float searchSpace, float beadradius);
    PointSetModel(float beadradius, float xaxis, float yaxis, float zaxis);
    PointSetModel(float beadradius, float height, float radius, std::string sym); //cylnder
    PointSetModel(float beadradius, float xaxis, float yaxis, float zaxis, std::string sym);
    PointSetModel(float searchSpace, float beadradius, std::string sym); //spherical
    PointSetModel(std::string helicalfile, float beadradius); // helical
    PointSetModel(float beadradius, float height, float radius); // cylinder no symmetry
    PointSetModel(std::string maskfile, float searchSpace, float beadradius); // masked search space
    PointSetModel(float beadradius, float xmin, float ymin, float zmin, float xmax, float ymax, float zmax, std::string sym);
    PointSetModel(std::vector<vector3> &centered_coordinates, float searchSpace, float beadradius);
    PointSetModel(std::string maskfile, float searchSpace, float beadradius, std::string sym, bool maskIt);

    Bead * getBead(int i) { return &beads[i];}

    void setSymmetryParameters(std::string sym);
    void createSubUnitAngles();
    float getZaxis(){ return zaxis;}
    float getXaxis(){ return xaxis;}
    float getYaxis(){ return yaxis;}

    float getNeighborCutOffLimit(){ return cutOffNeighbor;}
    float getBeadVolume() const {return bead_volume;}
    std::string getSymmetry(){ return symmetry;}

    unsigned int getNumberOfSubUnits() const;

    unsigned int getSymmetryIndex(){ return symmetryIndex;}
    std::string getSymmetryGroup(){ return symmetryGroup;}
    float * getSubUnitAnglesCos(){ return subUnitAnglesCos.data();}
    float * getSubUnitAnglesSin(){ return subUnitAnglesSin.data();}
    vector3 transformVectorBySymmetry(unsigned int subunitIndex, const vector3 & vec);
    void transformCoordinatesBySymmetry(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates);
    void transformCoordinatesBySymmetryPreCalc(const unsigned int subunitIndex, const unsigned int workingLimit, unsigned int &startCount, std::vector<vector3> &coordinates);

    void centerLatticeModel(unsigned int * workingLimit, std::vector<unsigned int> & indices, std::set<unsigned int> & hullPts);

    void createUniverse(bool useSphericalModel);
    void createDistancesAndConvertToSphericalCoordinates();
    void createUniverseFromMaskSym(std::string filename);
    void createUniverseFromMask(std::string filename);

    unsigned short int getMaxBin(PofRData * pData);
    float getBeadRadius(){ return bead_radius;}

    unsigned int getTotalDistances(){ return totalDistances;}
    unsigned int getTotalNumberOfBeadsInUniverse() const {return number_of_beads;}
    unsigned int getSizeOfNeighborhood(){ return sizeOfNeighborhood;}
    unsigned int getNeighborLimit() { return neighborLimit;}

    void setAverageNumberOfContactsInModel(float number){ this->averageNumberOfContactsInModel = number; }


    bool isUseDirectMethod() const { return useDirectMethod; }

    unsigned short int populateBins(PofRData * pData);
    unsigned short int * getPointerToBins(){ return &bins[0];}

    /*
     * this is sorted!
     */
    void copyStartingModelIntoVector(std::vector<unsigned int> & vec){std::copy(startingSetVector.begin(), startingSetVector.end(), vec.begin());}

    std::vector<unsigned int>::iterator getPointerToNeighborhood(unsigned int index);
    const unsigned int * getDirectPointerToNeighborhood(){ return this->neighbors.data(); }

    void pruneUniverseToAsymmetricUnit();

    std::string createHeader(float dkl, Anneal * annealedObject, PofRData *pData, unsigned int totalSteps, unsigned int workingNumber, float volume, float averageContacts);

    void writeModelToFile(const unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string & nameOf, unsigned int steps);
    void writeModelToFileFlipped(const unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string & nameOf, unsigned int steps);
    std::string writeBasicSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::string name);
    void writeSubModelToFile(unsigned int startIndex, unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::string nameOf);
    std::string writeModelToFile2(float dkl, unsigned int workingLimit, std::vector<unsigned int> &selectedBeads, std::vector<unsigned int> &pofrModel, std::string nameOf, Anneal *annealedObject, PofRData *pData, unsigned int steps, float volume, float averageContacts);
    void writeSetToFile(std::set<unsigned int> &selectedBeads, std::string & nameOf);
    std::string writeSymModelToFile(float dkl, unsigned int workingLimitS, std::vector<unsigned int> &beads, std::vector<unsigned int> &pofrModel, std::string name, Anneal *annealedObject, PofRData *pData, unsigned int totalSteps, float volume, float averageContacts);
    std::string writeModelToFileBare(float dkl, unsigned int workingNumber, std::vector<unsigned int> &selectedBeads, std::vector<unsigned int> &pofrModel, std::string nameOf, Anneal *annealedObject, unsigned int steps, float volume, float averageContacts);

    void setStartingSet(std::vector<unsigned int> & beads){
        total_in_working_universe = beads.size();
        startingSetVector.resize(total_in_working_universe);
        std::copy(beads.begin(), beads.end(), startingSetVector.begin());
    }

    /**
     * holds the number of beads, volume is detemined by multiplying by bead volume
     */
    void setBeadAverageAndStdev(float number_of_beads, float stdev) {
        beadAverage = number_of_beads;
        beadStDev = stdev;
    }

    void setStartingWorkingLimit(unsigned int limit){ startingWorkingLimit = limit;}
    unsigned int getStartingWorkingLimit(){return startingWorkingLimit;}

    unsigned int getBaseWorkingLimit(){ return baseWorkingLimit;}

    float getVolumeAverage(){return beadAverage;}
    float getVolumeStdev(){return beadStDev;}
    void setCVXHullVolume(float volume){ this->cvx_volume = volume; }
    float getCVXHullVolume () const { return this->cvx_volume; }


    unsigned int getTotalInSeed(){ return this->total_in_seed; }



    void createSeedFromPDB(std::string filename, PofRData * pData, unsigned int totalBins, std::vector<double> * pdbPr);


    void createSeedFromPDBSym(std::string filename, PofRData * pData, unsigned int totalBins, std::vector<double> * pdbPr);

    void copySeedVectorIntoSet(std::set<unsigned int> & pSet){ for(auto & val : seed_indices){ pSet.insert(val);};}

    const std::vector<unsigned int>::const_iterator getSeedBegin() const { return seed_indices.cbegin(); }
    const std::vector<unsigned int>::const_iterator getSeedEnd() const { return seed_indices.cend(); }
    const std::set<unsigned int>::const_iterator getReducedSeedBegin() const { return reduced_seed.cbegin(); }
    const std::set<unsigned int>::const_iterator getReducedSeedEnd() const { return reduced_seed.cend(); }

    void setReducedSeed(unsigned int limit, std::vector<unsigned int> &reduced){ //reduced_seed.resize(limit);
        total_in_reduced_seed = limit;

        reduced_seed.clear();
        for(unsigned int i=0; i<limit; i++){
            reduced_seed.insert(reduced[i]);
        }
    }

    inline void printBeadsFromSet(std::set<unsigned int> &beadIDs){

        std::string residue_index;
        for(auto it=beadIDs.begin(); it!=beadIDs.end(); ++it){
            int n = *it;
            residue_index = std::to_string(n+1);
            printf("%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, "CA", "ALA", "A", residue_index.c_str(), beads[n].getX(), beads[n].getY(), beads[n].getZ() );
        }
    }

    float getRadiusOfUniverse(){ return radius_of_universe;}

    void clearNeighborsAndBins(){
        bins.clear();
        neighbors.clear();
    }
};


#endif //IKETAMA_MODEL_H
