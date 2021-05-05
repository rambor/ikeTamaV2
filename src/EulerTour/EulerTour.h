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

class PointSetModel;

/**
 *   Class is used to monitor the graph connectivity
 *   EulerTour eulerTour(beginIt, subUnitWorkingLimit, pModel);
 *   currentNumberOfComponents = eulerTour.getNumberOfComponents();
 */

class EulerTour {

    struct find_Node_by_key : std::unary_function< Node *, bool>{
        find_Node_by_key(unsigned int keyToFind) : key(keyToFind){}
        bool operator () (Node * p) { return p->getKey() == key; }
    private:
        unsigned int key;
    };

    unsigned int totalComponents;
    std::map<unsigned int, Node> nodes; // use shared pointer? make instance on heap
    std::map<unsigned int, std::list< Node *> > tours; // key is the root of the tour

    void createInitialTour(unsigned int workingLimit, PointSetModel *pModel, std::vector<unsigned int>::const_iterator beginIt);
    bool addToTour(unsigned int nodeToAdd);

    void createSubTour(Node * pNode, std::list< Node * > * subTourToLoad);
    void rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);
    void rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);

    void printList(std::string text, std::list< Node * > * list);
    //void resetRootNodesInSubTour(std::list<Node *> * subTour);
    void resetRootNodesInSubTourOfTour(std::list<Node *> * subTour);
    bool validateTour(std::list<Node *> * tourtocheck);

public:
    // pass in iterator to beginning of selected lattice points that are sorted upto WorkingLimit
    // want to use this Class to determine connectivity of a graph
    // Connectivity is the numberOfComponents
    //
    EulerTour();
    EulerTour(std::vector<unsigned int>::const_iterator beginIt, const unsigned int workingLimit, PointSetModel *pModel);
    EulerTour(std::vector<unsigned int> & indices, const unsigned int workingLimit, PointSetModel *pModel);
    EulerTour(std::set<unsigned int> & beginIt, PointSetModel *pModel);

    EulerTour(const EulerTour &e2);

    ~EulerTour(){
        for(auto & tour : tours){
            for(auto & ele : tour.second){
                ele = nullptr;
            }
            tour.second.clear();
        }
        tours.clear();
        nodes.clear();
    }

    void swap(EulerTour & other) noexcept {
        using std::swap;
        std::swap(nodes, other.nodes);
        other.totalComponents = totalComponents;
        std::swap(tours, other.tours);
    }


    // copy assignment
    EulerTour & operator=(const EulerTour & tourtomove) {

        EulerTour tmp(tourtomove);
        tmp.swap(*this);

        return *this;
    }

    // terrible move constructor
    EulerTour(EulerTour && moveit) noexcept {
        totalComponents = moveit.totalComponents;

        for(auto & node : moveit.nodes){
            nodes.emplace(node.first, Node(node.second));
        }


        for(auto & e2Tour : moveit.tours){ // making a new copy

            auto pair = tours.emplace(e2Tour.first, std::list< Node *>()); // make new tour
            std::list<Node *> * pListIt = &(*pair.first).second;

            for(auto & lit : e2Tour.second){
                unsigned int keytocopy = lit->getKey();
                pListIt->push_back(&(*nodes.find(keytocopy)).second);
            }
        }

        nodes.clear();
        tours.clear();
    }

    unsigned int addNode(unsigned int latticePoint, PointSetModel *pModel);
    unsigned int removeNode(unsigned int indexOfNode);

    unsigned int getNumberOfComponents(){return totalComponents;}

    unsigned int newTour(std::vector<unsigned int>::iterator beginIt, unsigned int workingLimit, PointSetModel *pModel);

    void createInitialTour(unsigned int workingLimit, PointSetModel *pModel, std::vector<unsigned int> &indices);

    bool validateList(std::string text);
    bool validateNodes(std::string str);

    void printTourInfo(){
        for(auto & tour : tours){
            std::cout << "TOUR : " << tour.first << " SIZE : " << tour.second.size() << std::endl;
        }
    }

    const std::list<Node *> * const getPointerToLargestTour(){
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

    std::map<unsigned int, std::list< Node *> > * getTours(){ return &tours;} // key is the root of the tour
};



#endif //IKETAMA_EULERTOUR_H
