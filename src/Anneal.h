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

#ifndef IKETAMA_ANNEAL_H
#define IKETAMA_ANNEAL_H

#include "thirdparty/Eigen/Core"
#include "thirdparty/Eigen/Geometry"
#include "thirdparty/power_sasa.h"

#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <random>

#include "PointSetModel.h"
#include "KDE.h"
#include "../EulerTour/EulerTour.h"
#include "Component.h"
#include <sastools/utils.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <libqhull/qhull_a.h>
#ifdef __cplusplus
}
#endif

class Component; // forward declaration
struct ProbabilityBead;

class Anneal {
    struct Trial{
        double value;
        std::vector<unsigned int> indices;
        Trial() = default;

        Trial(double val, std::vector<unsigned int>  & vec) : value(val) {
            //indices = std::vector<unsigned int>(vec);
            indices = std::move(vec);
        }

        bool operator<(const Trial & a) const{
            return value < a.value;
        }
    };

//    struct find_pt_by_key : std::unary_function< KDE::ProbabilityBead, bool>{
//        find_pt_by_key(unsigned int keyToFind) : key(keyToFind){}
//        bool operator () (KDE::ProbabilityBead p) { return p.index == key; }
//    private:
//        unsigned int key;
//    };

    float percentAddRemove, beta, eta, lambda, alpha, mu, asaAcceptanceRate, complementASAAcceptanceRate, intASAAcceptanceRate, intComplementASAAcceptanceRate;
    float highT, highTempStartForCooling, interconnectivityCutOff;
    double targetContacts=1.99;
    unsigned int lowerV, upperV, highTempRounds, ccmultiple;
    bool isRefine=false;

    std::string filenameprefix;
    std::vector<double> contactsDistribution;

//    std::vector<Component> components;
    unsigned int totalComponents=0, distributionlimit=13;
    float contactCutOff, violation_limit, probe_radius, delta_r = 1.4f;
    unsigned int maxbin, totalBins;

    std::vector<Component> components;



public:
    Anneal(float highT,
           float percent,
           unsigned int highTempRounds,
           std::string prefix,
           float alpha,
           float beta,
           float eta,
           float lambda,
           float mu,
           unsigned int multiple,
           float accRate,
           float interconn
    );


    float getAsaAcceptanceRate(){ return asaAcceptanceRate;}

    void updateASAConstantTemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp);
    void updateASATemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp);


    void populateContactsDistribution(std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, PointSetModel *pModel);

    float surfaceToVolume(const unsigned int total, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates);
    float calculateCVXHULLVolume(char * flags, std::vector<unsigned int> *bead_indices, unsigned int upTo, PointSetModel *pModel);

    unsigned int numberOfContactsFromSet(std::set<unsigned int> *beads_in_use, PointSetModel *pModel, unsigned int selectedIndex);
    unsigned int numberOfContactsFromSetSym(std::set<unsigned int> *beads_in_use, PointSetModel *pModel, unsigned int selectedIndex);

    void fillPrBinsAndAssignTotalBin(PointSetModel * pModel, PofRData * pData);

    inline void beadToPoint(pointT *ptestPoint, Bead *pBead) {
        ptestPoint[0] = pBead->getX();
        ptestPoint[1] = pBead->getY();
        ptestPoint[2] = pBead->getZ();
    }

    unsigned int getMaxBin(){ return maxbin;}
    unsigned int gettotalBinsDerivedFromData(){ return totalBins;}
    unsigned int getTotalBins(){ return totalBins;}

    float getHighTempStartForCooling(){return highTempStartForCooling;}
    float getPercentAddRemove(){ return percentAddRemove;}
    float getAlpha(){ return alpha;}
    float getLambda(){ return lambda;}
    float getMu(){ return mu;}
    float getEta(){ return eta;}
    float getBeta(){ return beta;}


    void getHullPoints(std::set<unsigned int> &hullpts, std::set<unsigned int> &beads_in_use, PointSetModel * pModel);

    void ceMapOptimization(PointSetModel *pModel, PofRData *pData, unsigned int topN, unsigned int minN, unsigned int maxN, std::vector<ProbabilityBead> & lattice);
    void ceMapOptimizationSym(PointSetModel *pModel, PofRData *pData, unsigned int topN, const unsigned int minN, const unsigned int maxN, std::vector<ProbabilityBead> & lattice);
    void calculateModelPrDistributionDirect(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, const unsigned int workingLimit, PointSetModel *pModel, PofRData *pData);
    void calculateModelPrDistributionSymCE(std::vector<unsigned int> *subUnit_indices,
                                           std::set<unsigned int> & subUnit_indices_tree,
                                           std::map<unsigned int, std::vector<vector3> > &map,
                                           std::vector<double> & contactsDistributionOfModel,
                                           std::vector<unsigned int> *binCount,
                                           const unsigned int indicesWorkingLimit,
                                           unsigned int &violations,
                                           PointSetModel *pModel, PofRData *pData);


    bool createInitialModelCVXHull(PointSetModel *pModel, PofRData *pData, std::string name);
    std::string refineHomogenousBodyASAHybridEx(PointSetModel *pModel, PofRData *pData, std::string outputname);

    bool createInitialModelSymmetry(PointSetModel *pModel, PofRData *pData);
    std::string refineSymModel(PointSetModel *pModel, PofRData *pData, std::string nameTo);

    bool initializeModelToRefineSym(PointSetModel *pModel, PofRData *pData, std::string name, std::string PDBFilename);
    std::string refineSymModelRefine(PointSetModel *pModel, PofRData *pData, std::string nameTo);

    bool createInitialModelCVXHullDirect(PointSetModel *pModel, PofRData *pData, std::string name);
    std::string refineHomogenousBodyASAHybridDirect(PointSetModel *pModel, PofRData *pData, std::string outputname);

    bool createSeedFromPDB(PointSetModel *pModel, PofRData *pData, std::string name, std::string PDBFilename);
    std::string refineHomogenousBodyASACVXSeeded(PointSetModel *pModel, PofRData *pData, std::string outputname);

    void calculateModelParametersSymmetry(std::set<unsigned int> *subUnit_indices, PointSetModel *pModel);

    void recalculateDeadLimit(unsigned int workingLimit, std::vector<unsigned int> &bead_indices, std::set<unsigned int> &hull, PointSetModel * pModel);
    void calculateModelPrDistribution(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, unsigned int workingLimit, unsigned int totalBeadsInSphere, PointSetModel *pModel);
    void calculateModelPrDistributionSym(std::vector<unsigned int> *subUnit_indices, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, unsigned int &violation,  PointSetModel *pModel, PofRData *pData);

    void populateLayeredDeadlimitUsingSet(std::set<unsigned int> & beads_in_use, std::set<unsigned int> & hull, PointSetModel * pModel);

    void addToContactsDistribution(unsigned int addMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, PointSetModel *pModel);
    double connectivityPotential(unsigned int numberOfComponents);

    void printParameters(std::vector<float> * accept, std::vector<double> * temp, std::vector<float> * divergence, std::vector<unsigned int> * wl);

    void removeFromContactsDistribution(
            unsigned int removeMe,
            std::vector<double> & distribution,
            std::set<unsigned int> *beads_in_use,
            PointSetModel *pModel);

    void removeLatticePositionToModelDirect(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned int * pWorkingLimit,
            const unsigned int * pLatticePointToRemove,
            PointSetModel * pModel,
            PofRData * pData);

    unsigned int getConnectedNeighborFromSet(
            std::set<unsigned int> *beads_in_use,
            PointSetModel *pModel,
            unsigned int & selectedIndex);

    unsigned int getUseableNeighborFromSet(std::set<unsigned int> *beads_in_use,
                                                   PointSetModel *pModel,
                                                   unsigned int & selectedIndex);

    unsigned int getUseableNeighborFromSetCE(std::set<unsigned int> *beads_in_use,
                                           PointSetModel *pModel,
                                           unsigned int & selectedIndex);


    unsigned int getUseableNeighborWeighted(std::set<unsigned int> *beads_in_use,
                                                    PointSetModel *pModel,
                                                    unsigned int & selectedIndex);

    float getProbabilityOfSelectedIndex(std::set<unsigned int> *beads_in_use,
                                                       PointSetModel *pModel,
                                                       unsigned int & selectedIndex);

    unsigned int removeFromPrSym(
            unsigned int const &removeMeSubUnitIndex,
            std::vector<unsigned int> & beadsInUse,
            unsigned int const &workingLimit,
            std::vector<unsigned int> & prBins,
            PointSetModel *pModel,
            PofRData *pData);


    unsigned int addToPrSym(unsigned int addMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, PointSetModel *pModel, PofRData *pData);

    unsigned int getViolations(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, PointSetModel *pModel);
    unsigned int getViolationsFromSet(std::set<unsigned int> *subUnit_indices, std::vector<double> & contactsDistributionOfModel, unsigned int indicesWorkingLimit, PointSetModel *pModel);

    void addLatticPositionToModelDirect(std::vector<unsigned int> * pIndices,
                                        unsigned int * pWorkingLimit,
                                        std::vector<unsigned int>::iterator * pItIndex);

    void addToPrDirect(const unsigned int addMe,
                       std::vector<unsigned int> & beadsInUse,
                       const unsigned int upperLimit,
                       std::vector<unsigned int> & prBins,
                       PointSetModel * pModel,
                       PofRData * pData );

    void restoreAddingFromBackUpDirect(
            unsigned int * pWorkingLimit,
            std::vector<unsigned int> * pBinCountBackUp,
            std::vector<unsigned int> * pBinCount);

    void restoreRemovingLatticePointFromBackUpDirect(
            unsigned int * pWorkingLimit,
            std::vector<unsigned int> * pBinCountBackUp,
            std::vector<unsigned int> * pBinCount);

    void removeFromPrDirect(unsigned int removeMe,
                            std::vector<unsigned int> & beadsInUse,
                            const unsigned int upperLimit,
                            std::vector<unsigned int> & prBins,
                            PointSetModel * pModel,
                            PofRData * pData);

    void fillPrBinsAndAssignTotalBinDirect(PointSetModel * pModel, PofRData * pData);

    float calculateKLDivergenceAgainstPDBPR(std::vector<unsigned int> &modelPR, std::vector<double> &targetPR);

    void setContactsDistribution(std::vector<unsigned int> & contactsDistributionSeed);
    bool setAnchorPoints(std::string anchorFileName, std::string pdbFile, PointSetModel *pModel);

    bool canRemoveIfNotCenteredAnchor(unsigned int index);

    float connectivityPotentialPhases(unsigned int mainConnectivity);

    double estimateMagnitudeOfDifferenceContactsPotential(unsigned int workingLimit,
                                                          std::vector<unsigned int> & bead_indices,
                                                          std::vector<unsigned int> & binCount,
                                                          PointSetModel * pModel);
    /**
    *
    */
    inline void addLatticPositionToModel(std::vector<unsigned int> * pIndices,
                                         unsigned int * pWorkingLimit,
                                         std::vector<unsigned int>::iterator * pItIndex){


        //std::copy(pIndices->begin(), pIndices->end(), pBackUpState->begin()); // make backup
        // make the swap at the border (workingLimit)
        std::iter_swap(pIndices->begin() + *pWorkingLimit, *pItIndex); // this swaps to working position and changes backup
        // increment workingLimit to include new position
        *pWorkingLimit += 1;
        std::sort(pIndices->begin(), pIndices->begin() + *pWorkingLimit);
    }

    /**
     * addMe must be present in beadsInUse
     * pBin is a pointer to the vector of binned distances in Model
     * prBin is the vector holding counts per bin
     * upperLimit is usually the working Limit of the set
     */
    inline void addToPr(const unsigned int addMe, const std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, const unsigned short int * const pBin, const unsigned int & totalBeadsInUniverse, std::vector<unsigned int> & prBins ){

        // beadsInUse must be sorted
        unsigned int row;
        unsigned int row2;


        unsigned int i=0;
        // add down column
        const unsigned int * const ptr = beadsInUse.data();

        while(ptr[i] < addMe && i < upperLimit){
            row = ptr[i];
            row2 = (row*totalBeadsInUniverse - (row*(row+1)/2)) + (addMe - row);// to prevent underFlow subtract 1 last
            prBins[ *(pBin + row2 - 1 ) ]++;
            i++;
        }
        // out of first loop, i = addMe
        i++;
        // Add across row
        row2 = addMe*(unsigned int)totalBeadsInUniverse - addMe*(addMe+1)/2 - addMe;
        while (i < upperLimit){
            prBins[ *(pBin + (row2 + ptr[i]) - 1 ) ]++;
            i++;
        }
    }

    inline void restoreAddingFromBackUp(std::vector<unsigned int> * pIndices,
                                        const std::vector<unsigned int> * pBackUpState,
                                        unsigned int * pWorkingLimit,
                                        const std::vector<unsigned int> * pBinCountBackUp,
                                        std::vector<unsigned int>::iterator * pBinCountBegin){

        std::copy(pBackUpState->begin(), pBackUpState->end(), pIndices->begin());
        *pWorkingLimit -= 1;
        std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), *pBinCountBegin); //copy to bin count
    }

    /**
     *
     */
    inline void removeLatticePositionToModelByIndex(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned short int * const pBin,
            unsigned int * pWorkingLimit,
            unsigned int totalBeadsInSphere,
            const unsigned int latticeIndexToRemove){

        const auto itIndex = bead_indices.begin() + latticeIndexToRemove;
        // remove original from P(r)
        // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
        removeFromPr(*itIndex, bead_indices, *pWorkingLimit, pBin, totalBeadsInSphere, pBinCount);
        // reduce the workingLimit
        // if wl = 10
        // 0 1 2 3 4 5 6 7 8 9 10
        // remove 4, wl-=1
        // 0 1 2 3 9 5 6 7 8 4 10
        // sort to 9
        // 0 1 2 3 5 6 7 8 9 4 10
        //
        *pWorkingLimit -= 1;
        // swap selected point to the workingLimit
        std::iter_swap(itIndex, bead_indices.begin() + *pWorkingLimit);
        // still need to sort, swap changes the order
        // to save some time, sort should only be from swapped point to workingLimit
        std::sort(bead_indices.begin(), bead_indices.begin() + *pWorkingLimit);
    }


    /**
      * requires a sorted beadsInuse vector
      * removeMe is the index of the bead from origination
    */
    inline void removeFromPr(const unsigned int removeMe, std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, unsigned short int * const pBin, const unsigned int & totalBeadsInUniverse, std::vector<unsigned int> & prBins ){
        // beads in Use is sorted
        unsigned int row;
        unsigned int row2;

        unsigned int i=0;

        //const int * ptr = &(beadsInUse.front());
        const unsigned int * ptr = beadsInUse.data();

        while(ptr[i] < removeMe && i < upperLimit){
            row = ptr[i];
            row2 = (row*totalBeadsInUniverse - (row*(row+1)/2)) + (removeMe - row);// - 1;
            prBins[ *(pBin + row2 - 1 ) ]--;
            i++;
        }

        // when loop ends beadsInUse[i] => removeMe
        i++;
        // remove across row
        row2 = removeMe*totalBeadsInUniverse - removeMe*(removeMe+1)/2 - removeMe; // assume this is always positive
        while (i < upperLimit){
            //prBins[ *(pBin + row2 + beadsInUse[i] - 1 ) ]--;
            prBins[ *(pBin + (row2 + ptr[i]) - 1 ) ]--;
            i++;
        }
    }

    inline void restoreRemovingLatticePointFromBackUp(std::vector<unsigned int>::iterator * pBeginIt,
                                                      unsigned int * pWorkingLimit,
                                                      std::vector<unsigned int> * pBinCountBackUp,
                                                      std::vector<unsigned int>::iterator * pBinCountBegin){
        // value we want to recover is at wl = 9
        // wl += 1 => 10
        // 0 1 2 3 9 5 6 7 8 4 10
        // sort to wl
        // 0 1 2 3 4 5 6 7 8 9 10
        //
        *pWorkingLimit += 1;
        // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
        // sorting is n*log(n) with n = 200?  should be much smaller
        std::sort(*pBeginIt, *pBeginIt+*pWorkingLimit); // swapped index is at workingLimit
        std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), *pBinCountBegin); //copy to bin count
    }


    /**
   *
   * calculate Kulback Liebler divergence of contacts distribution
   * @param distribution
   * @return
   */
    inline double calculateAvgContactsPotential(std::vector<unsigned int> &binCount, unsigned int wl) {

        double avgCon = 2.0d*binCount[0]/(float)wl;
        double baseline = 0.0000001;

        double diff = 0;
//        diff = (avgCon-targetContacts);
//        diff *= diff;

//        double bfactor = -0.3719;
//        diff = 1.0 - std::exp(-bfactor*(avgCon-targetContacts));
//        diff *= diff;
//        diff *= 0.00001;

        //diff += baseline;

//        if (avgCon < targetContacts){ // lower bound, can't be less than
//            diff = (avgCon-targetContacts);
//            diff *= diff;
//            diff += baseline;
//        } else {
//            double slope = (0.0000001d - 1.0d)*baseline/(12.0 - targetContacts); // drop tenth of baseline over 12 units
//            double intercept = baseline - slope*targetContacts;
//            diff = slope*avgCon + intercept;
//        }


//        if (avgCon < targetContacts){ // lower bound, can't be less than
//            diff = (avgCon-targetContacts);
//            diff *= diff;
//            diff += baseline;
//        } else {
//            double slope = 0.001d*baseline/(12.0 - targetContacts); // drop tenth of baseline over 12 units
//            double intercept = baseline - slope*targetContacts;
//            diff = slope*avgCon + intercept;
//        }



//        double upperlimit = 5.1;
//        if (avgCon < targetContacts){ // lower bound, can't be less than
//            diff = (avgCon-targetContacts);
//            diff *= diff;
//            diff += baseline;
//        } else if (avgCon > upperlimit){
//            diff = (avgCon-upperlimit);
//            diff *= diff;
//            diff += baseline;
//        } else {
//            diff += baseline;
//        }


//        if (avgCon > targetContacts){ // upper bound on contacts, can't be any greater than target
//            diff = (avgCon-targetContacts);
//            diff *= diff;
//            diff += baseline;
//        } else {
//            double delta = 0.0001;
//            double slope = delta*baseline/targetContacts; // increase by tenth of baseline
//            diff = slope*avgCon + (1.0 - delta)*baseline;
//        }

        return diff;
    }


    /**
       *
       * calculate Kulback Liebler divergence of contacts distribution
       * @param distribution
       * @return
       */
    inline double calculateKLDivergenceContactsDistribution(std::vector<double> &distribution) {

        double kl=0.0d, totalCounts = 0.0, avg = 0.0;
        for(unsigned int i=0; i< 13; i++){
            totalCounts += distribution[i];
            avg += i*distribution[i];
        }

        double avgCon = avg/totalCounts;
        double diff = 0, baseline = 0.000001;

        if (avgCon < targetContacts){ // lower bound, can't be less than
            diff = (avgCon-targetContacts);
            diff *= diff;
            diff += baseline;
        } else {
            double slope = -0.001*baseline/(12.0 - targetContacts); // drop tenth of baseline over 12 units
            double intercept = baseline - slope*targetContacts;
            diff = slope*avgCon + intercept;
        }

//        for (unsigned int i=0; i < distributionlimit; i++){
//            // i know every value in working_probability up to zeroBin is nonzero
//            double prob = contactsDistribution[i];  // bounded by experimental Shannon Number
//            double tempvalue = distribution[i];
//            if (prob > 0){
//                if (tempvalue > 0){
//                    kl += prob * std::log(prob / tempvalue * totalCounts);
//                } else { // if tempvalue is Zero (empty)
//                    kl += 1.1;
//                }
//            }
//        }

//        double percent = distribution[1]/totalCounts;
        double diff2 = 0;//distribution[1]/totalCounts - contactsDistribution[1]; // number of beads that are single connected

        for(unsigned int i=5; i< 13; i++){
            diff2 += distribution[i]/totalCounts;
        }

//        return diff*diff;
//         can't be any more than
//        if (percent > contactsDistribution[1]){
//            double difft = percent - contactsDistribution[1];
//            diff2 = difft*difft + 1.1*baseline;
//        } else {
//            double slope = 0.1*baseline/contactsDistribution[1]; // increase by tenth of baseline
//            diff2 = (slope*percent + baseline);
//        }

        //return kl;
        //return morse*morse;
        return diff2*diff2;
        //return diff;
    }

    /**
 *
 */
    inline void removeLatticePositionToModel(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned short int * const pBin,
            unsigned int * pWorkingLimit,
            unsigned int totalBeadsInSphere,
            const unsigned int * pLatticePointToRemove){

        auto pBeginIt = bead_indices.begin();
        auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
        // remove original from P(r)
        // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
        removeFromPr(*pLatticePointToRemove, bead_indices, *pWorkingLimit, pBin, totalBeadsInSphere, pBinCount);
        // reduce the workingLimit
        // if wl = 10
        // 0 1 2 3 4 5 6 7 8 9 10
        // remove 4, wl-=1
        // 0 1 2 3 9 5 6 7 8 4 10
        // sort to 9
        // 0 1 2 3 5 6 7 8 9 4 10
        //
        *pWorkingLimit -= 1;
        // swap selected point to the workingLimit
        std::iter_swap(itIndex, pBeginIt + *pWorkingLimit);
        // still need to sort, swap changes the order
        // to save some time, sort should only be from swapped point to workingLimit
        std::sort(pBeginIt, pBeginIt + *pWorkingLimit);
    }

    /**
    *
    *
    */
    inline unsigned int removeLatticePositionToModelSym(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & modelPrBins,
            unsigned int * pWorkingLimit,
            const unsigned int * pLatticePointToRemove,
            PointSetModel * pModel,
            PofRData *pData){
        auto pBeginIt = bead_indices.begin();
        auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
        // remove original from P(r)
        // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
        unsigned int violations = removeFromPrSym(*pLatticePointToRemove, bead_indices, *pWorkingLimit, modelPrBins, pModel, pData);
        // reduce the workingLimit
        // if wl = 10
        // 0 1 2 3 4 5 6 7 8 9 10
        // remove 4
        // 0 1 2 3 9 5 6 7 8 4 10
        //
        // 0 1 2 3 5 6 7 8 9 4 10
        //
        *pWorkingLimit -= 1;
        // swap selected point to the workingLimit
        std::iter_swap(itIndex, pBeginIt + *pWorkingLimit);
        // still need to sort, swap changes the order
        // to save some time, sort should only be from swapped point to workingLimit
        //std::sort(pBeginIt, pBeginIt + *pWorkingLimit);
        return violations;
    }

    // coupon collector's problem
    inline unsigned int getCoupons(unsigned int numberOfPoints){
        auto couponbase = (double) numberOfPoints;
        return (unsigned int)(couponbase*std::log(couponbase) + 0.5772156649*couponbase + 0.5d);
    }



};


#endif //IKETAMA_ANNEAL_H
