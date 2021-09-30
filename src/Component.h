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

#ifndef IKETAMA_COMPONENT_H
#define IKETAMA_COMPONENT_H

#include <string>
#include <vector>
#include <cmath>
#include <regex>
#include <iostream>
#include <set>
#include "../EulerTour/EulerTour.h"
#include <sastools/Bead.h>
#include "PointSetModel.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <libqhull/qhull_a.h>
#ifdef __cplusplus
}
#endif


class Component {

private:

    std::string id;
    std::set<unsigned int> anchors{};            // can be empty
    std::set<unsigned int> centeredAnchors{};            // can be empty, anchor is bead that is centered on the residued

    std::vector<unsigned int> resids{};          // can be empty
    std::vector<std::string> chains{};  // can be empty
    std::vector<float> potentialFunction{};  // can be empty

    std::set<unsigned int> beads_in_use{};       // never empty
    std::set<unsigned int> best{};       // never empty

    std::vector<unsigned int> cvxPoints{};       // never empty
    std::vector<double> contactsDistributionOfModel{};
    EulerTour tour;

    std::set<unsigned int> hull{};

    unsigned int totalComponents;
    unsigned int targetNumberOfBeads=0;
    float targetVolume, invTargetVolume;
    PointSetModel * pModel = nullptr;  // for the euler tour
    float currentCVXVolume;
    bool empty = true;
    double totalContactSum;
    unsigned int anchorCount=0;
    float percentageStep;

    double contactsKLDivergence;

public:

    Component(std::string id, float volume, PointSetModel *pModel);
    Component(std::string id, float volume, PointSetModel *pModel, bool contiguous);

    ~Component(){
        //delete pModel; // object is not handling memory of pModel
        pModel = nullptr;
    }

    Component(const Component &toCopy){
        id = toCopy.id;
        anchors = std::set<unsigned int>(toCopy.anchors);
        centeredAnchors = std::set<unsigned int>(toCopy.centeredAnchors);

        resids = std::vector<unsigned int>(toCopy.resids);
        chains = std::vector<std::string>(toCopy.chains);
        potentialFunction = std::vector<float>(toCopy.potentialFunction);

        beads_in_use = std::set<unsigned int>(toCopy.beads_in_use);
        best = std::set<unsigned int>(toCopy.best);

        cvxPoints = std::vector<unsigned int>(toCopy.cvxPoints);
        contactsDistributionOfModel = std::vector<double>(toCopy.contactsDistributionOfModel);

        tour = toCopy.tour;

        hull = std::set<unsigned int>(toCopy.hull);
        totalComponents = toCopy.totalComponents;
        targetNumberOfBeads = toCopy.targetNumberOfBeads;
        targetVolume= toCopy.targetVolume;
        invTargetVolume = toCopy.invTargetVolume;

        pModel = &*(toCopy.pModel); // should point to the same model

        currentCVXVolume = toCopy.currentCVXVolume;
        empty = toCopy.empty;

        totalContactSum = toCopy.totalContactSum;
        anchorCount = toCopy.anchorCount;

        percentageStep = toCopy.percentageStep;
        contactsKLDivergence = toCopy.contactsKLDivergence;
    }

    // copy assignment operator
    Component & operator=(const Component & dataToCopy) {
        if (this == &dataToCopy)
            return *this;

        Component tmp(dataToCopy); // make a copy
        tmp.swap(*this);
        return *this;
    }

    /**
     * Rule of 3.5, define copy, destructor and assignment operator
     * @param other
     */
    void swap(Component & other) {

        std::swap(pModel, other.pModel);

        other.anchors = std::move(anchors);
        other.centeredAnchors = std::move(centeredAnchors);

        other.resids = std::move(resids);
        other.chains = std::move(chains);
        other.potentialFunction = std::move(potentialFunction);
        
        other.beads_in_use = std::move(beads_in_use);
        other.best = std::move(best);
        
        other.cvxPoints = std::move(cvxPoints);
        other.contactsDistributionOfModel = std::move(contactsDistributionOfModel);

        other.tour = std::move(other.tour); // uses move operator in EulerTour
        other.hull = std::move(other.hull);
        
        other.totalComponents = other.totalComponents;
        other.targetNumberOfBeads = other.targetNumberOfBeads;
        other.targetVolume= other.targetVolume;
        other.invTargetVolume = other.invTargetVolume;
        other.currentCVXVolume = other.currentCVXVolume;
        other.empty = other.empty;

        other.totalContactSum = other.totalContactSum;
        other.anchorCount = other.anchorCount;

        other.percentageStep = other.percentageStep;
        other.contactsKLDivergence = other.contactsKLDivergence;
    }

    bool checkID(std::string str);
    float numberOfLatticePointsPotential(float value);

    unsigned int removeLatticePoint(unsigned int index);
    unsigned int addLatticePoint(unsigned int index); // add to beads_in_use and

    bool indexInUse(unsigned int index);

    void addResid(unsigned int index, std::string chain);

    std::string const getID() const { return id;}
    unsigned int getTotalResids() { return resids.size();}
    unsigned int getResidByIndex(unsigned int index) { return resids[index];}
    std::string getChainsByIndex(unsigned int index) { return chains[index];}
    void addAnchor(unsigned int index);

    // desired number of beads but exclude the anchors
    void setTargetNumberOfLatticePoints(float value);
    void addCenteredAnchors(unsigned int index);
    unsigned int getTargetNumberOfLatticePoints();

    bool anchorsEmpty(){ return empty;}

    std::set<unsigned int> * getAnchors(){ return &anchors; }
    std::set<unsigned int> * getCenteredAnchors(){ return &centeredAnchors; }
    std::set<unsigned int> * getBeadsInUse() { return &beads_in_use; }
    std::set<unsigned int> * getHull(){ return &hull;}
    std::vector<double> * getContactsDistribution(){ return &contactsDistributionOfModel;}

    unsigned int getTotalNumberOfBeads() { return beads_in_use.size(); }

    unsigned int getTotalNumberOfComponents(){ return tour.getNumberOfComponents();}

    bool inUse(unsigned int index){ return !(beads_in_use.find(index) == beads_in_use.end()); }
    bool isAnchor(unsigned int index);
    bool isCenteredAnchor(unsigned int index);
    float potential();
    void printAnchors();
    void printCenteredAnchors();
    void printContactsPerCenteredAnchor();

    void printSet();
    void printBest();
    void writeToFile(std::string nameOf);
    void writeAnchorsToFile(std::string nameOf);
    void writeCenteredAnchorsToFile();
    void populatePotential(float percentage);

    void copyBestToInUse();
    float calculateCVXVolume();
    void setBest();

    void setTotalContactSum(double value){ totalContactSum = value;}
    double getTotalContactSum(){ return totalContactSum;}
    void printConstraints();
    unsigned int getAnchorCount(){return anchorCount;}
    void setCurrentCVXVolume();
    void setCurrentCVXVolume(float value){ this->currentCVXVolume=value;}
    float getCurrentCVXVolume(){ return currentCVXVolume;}

    //float getCurrentVolumeEnergy(){ return std::abs(currentCVXVolume-targetVolume)*invTargetVolume;}
    float calculateVolumeEnergy(float testVolume){ return std::abs(testVolume-targetVolume)*invTargetVolume;}


    double setContactsKLDivergence(double val){ return contactsKLDivergence = val;}
    double getContactsKLDivergence(){ return contactsKLDivergence;}

};


#endif //IKETAMA_COMPONENT_H
