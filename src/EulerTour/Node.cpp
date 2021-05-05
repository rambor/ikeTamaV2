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

#include "Node.h"

/**
 * CAN NOT HAVE A NODE WITH NO KEY!
 */

Node::Node(unsigned int key) : key(key), rootNode(key){
    totalNeighbors = 0;
}


void Node::addNeighbor(Node * pNode) {

    adjacencyListIndex.insert(pNode->getKey()); // set can never have duplicated values
    totalNeighbors = (unsigned int)adjacencyListIndex.size();

    if (!(pNode->isNeighborPresent(this->getKey()))){
        pNode->addNeighbor(this);
    }
}



void Node::removeNeighbor(Node * pNode) {

    adjacencyListIndex.erase(pNode->getKey());

//    if (adjacencyListIndex.find(pNode->getKey()) != adjacencyListIndex.end()){
//        if (totalNeighbors == 1){
//            adjacencyList[0] = 0;
//            adjacencyListIndex.erase(pNode->getKey()); // remove from set
//            totalNeighbors=0;
//        } else {
//
//            for(unsigned int i=0;i < totalNeighbors; i++){
//                //std::cout << this->getKey() << " - " << pNode->getKey() << " " << adjacencyList[i]->getKey() << " " << (adjacencyList[i] == pNode) << std::endl;
//                if (adjacencyList[i] == pNode->getKey()){ // do they have the same address (comparing pointers)
//
//                    std::iter_swap(adjacencyList.begin()+i, adjacencyList.begin()+totalNeighbors - 1);
//
//                    adjacencyListIndex.erase(pNode->getKey()); // remove from set
//                    totalNeighbors--;
//                    break;
//                } else {
//                    std::cout << "Parent NODE : " << this->getKey() << " not found neighbor during remove => " << pNode->getKey() << std::endl;
//                }
//            }
//        }
//    }

    totalNeighbors = (unsigned int)adjacencyListIndex.size();
    // remove itself from pNode's neighborhood list
    if (pNode->isNeighborPresent(this->getKey())){
        pNode->removeNeighbor(this);
    }
}

//void Node::setFirst(std::list<Node>::iterator *pFirst) {
//    this->pFirst = pFirst;
//}
//
//void Node::setLast(std::list<Node>::iterator *pLast) {
//    this->pLast = pLast;
//}

bool Node::validate(){

    if (totalNeighbors != adjacencyListIndex.size()){
        std::cout << "Possible repeats "<< std::endl;
        std::cout << "ADJACENCYLISTINDEX (SET) : " << adjacencyListIndex.size() << " != " << totalNeighbors << std::endl;

        return false;
    }
    //this->printNeighbors();
    return true;
}


void Node::printNeighbors(){
//    int count=1;
//    std::cout << " SIZE OF ADJACENCY LIST : " << adjacencyList.size() << " | TOTAL NEIGHBORS : " << totalNeighbors << " "  << std::endl;
//
//    for(int i=0; i<totalNeighbors; i++){
//        std::cout << "      " << count << " NEIGHBOR OF " << key << " " << adjacencyList[i] << " " << totalNeighbors << std::endl;
//        count++;
//    }
    // print adjacency list
    this->printAdjacencyList("FROM PRINTNEIGHBORS");
}


void Node::printAdjacencyList(std::string text){
    std::cout << "      ADJACENCY LIST " << std::endl;
    int count=1;
//    for(auto it= adjacencyListIndex.begin(); it!=adjacencyListIndex.end(); ++it){
//        std::cout << "      " << count << " ADJACENCY LIST : " << key << " " << *it << std::endl;
//        count++;
//    }
    for(auto it : adjacencyListIndex){
        std::cout << this->getKey() << "   |   " << count << " ADJACENCY LIST : " << key << " " << it << std::endl;
        count++;
    }
}


// setting pointer to tour should override previous root node
// void Node::setPointerToTour(std::list < Node * > * pointer){ // point to a location in map
//    pointerToTour = pointer;
//    rootNode = (*pointerToTour).front()->getKey();
////    std::cout << key << "          ROOT : SETTING POINT TO TOUR ( SIZE : " << pointer->size() << " ) "<< rootNode << std::endl;  // no size means list is empty
//}

bool Node::isNeighborPresent(unsigned int index) {
    return (adjacencyListIndex.find(index) != adjacencyListIndex.end());
}

