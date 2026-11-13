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

#ifndef IKETAMA_EULERTOUR_H
#define IKETAMA_EULERTOUR_H


#include <string>
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include "Node.h"
#include <ctime>
#include <set>
#include <utility>
#include <iostream>

class PointSetModel;

/**
 *   Tracks the number of connected components (graph connectivity) of the selected lattice points.
 *   EulerTour eulerTour(beginIt, subUnitWorkingLimit, pModel);
 *   currentNumberOfComponents = eulerTour.getNumberOfComponents();
 *
 *   Same interface and same results as the Euler-tour implementation, different internals:
 *   every occupied lattice point carries a component label.
 *     addNode     merges the components of the new point's neighbours (smaller ones are relabelled)
 *     removeNode  if the removed point's neighbours are still connected through each other (checked locally,
 *                 within the neighbourhood), nothing changes; otherwise a breadth-first search runs from each
 *                 piece in turn and stops as soon as the pieces meet or all but one are exhausted, so the cost is
 *                 set by the smaller pieces, not the model size.
 *   The previous version rebuilt and searched std::list tours with std::map/std::set look-ups on every removal,
 *   which was the largest single cost in the annealing loops.
 *
 *   getTours()/getPointerToLargestTour() are kept for compatibility: they build, on demand, one list of Node*
 *   per component (each node listed once, not an Euler tour).
 */
class EulerTour {

    PointSetModel * pModel = nullptr;
    const unsigned int * neighborhood = nullptr;
    unsigned int sizeOfNeighborhood = 12, neighborLimit = 0, universe = 0;

    unsigned int totalComponents = 0;
    std::vector<int> label;                  // component label per lattice point, -1 if not in the set
    std::vector<unsigned int> compSize;      // size per label
    std::vector<int> freeLabels;
    unsigned int totalNodes = 0;

    // work space for the searches (epoch-stamped so nothing has to be cleared)
    std::vector<unsigned int> stamp;
    std::vector<unsigned short> owner;
    unsigned int epoch = 0;

    // compatibility views, built on demand
    std::map<unsigned int, Node> nodes;
    std::map<unsigned int, std::list< Node *> > tours;

    void attach(PointSetModel * model);
    int newLabel();
    void relabelComponent(unsigned int start, int from, int to);
    void labelAll(const std::vector<unsigned int> & points);
    void buildCompatibilityViews();
    unsigned int nextEpoch();

public:
    EulerTour();
    EulerTour(std::vector<unsigned int>::const_iterator beginIt, const unsigned int workingLimit, PointSetModel *pModel);
    EulerTour(std::vector<unsigned int> & indices, const unsigned int workingLimit, PointSetModel *pModel);
    EulerTour(std::set<unsigned int> & beginIt, PointSetModel *pModel);

    EulerTour(const EulerTour &e2) = default;
    EulerTour(EulerTour && moveit) noexcept = default;
    EulerTour & operator=(const EulerTour & tourtomove) = default;
    EulerTour & operator=(EulerTour && tourtomove) noexcept = default;
    ~EulerTour() = default;

    void swap(EulerTour & other) noexcept { std::swap(*this, other); }

    unsigned int addNode(unsigned int latticePoint, PointSetModel *pModel);
    unsigned int removeNode(unsigned int indexOfNode);

    unsigned int getNumberOfComponents(){return totalComponents;}

    unsigned int newTour(std::vector<unsigned int>::iterator beginIt, unsigned int workingLimit, PointSetModel *pModel);

    void createInitialTour(unsigned int workingLimit, PointSetModel *pModel, std::vector<unsigned int> &indices);

    // consistency checks against a full recount: validateList returns true on ERROR, validateNodes true if OK (as before)
    bool validateList(std::string text);
    bool validateNodes(std::string str);
    unsigned int countComponentsBruteForce() const;

    void printTourInfo(){
        buildCompatibilityViews();
        for(auto & tour : tours){
            std::cout << "TOUR : " << tour.first << " SIZE : " << tour.second.size() << std::endl;
        }
    }

    const std::list<Node *> * const getPointerToLargestTour(){
        buildCompatibilityViews();
        unsigned int max = 0;
        std::list< Node *> * pList = nullptr;
        for(auto & tour : tours){
            if (tour.second.size() > max){
                max = tour.second.size();
                pList = &tour.second;
            }
        }
        return pList;
    }

    std::map<unsigned int, std::list< Node *> > * getTours(){ buildCompatibilityViews(); return &tours;}
};



#endif //IKETAMA_EULERTOUR_H
