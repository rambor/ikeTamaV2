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

#ifndef IKETAMA_NODE_H
#define IKETAMA_NODE_H

#include <string>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <cstdlib>
#include <iostream>

class Node {

private:
    unsigned int key{};
    unsigned int totalNeighbors = 0;
    unsigned int rootNode; // single node on its own is itself

    std::set< unsigned int > adjacencyListIndex;
    bool accessed = false;
    void printAdjacencyList(std::string text);

public:
    //Node();
    explicit Node(unsigned int key);

    // copy constructor
    Node(const Node &n2){
        this->key = n2.key;
        this->rootNode = n2.rootNode;
        for(auto & it : n2.adjacencyListIndex){
            this->adjacencyListIndex.insert(it);
        }
        this->totalNeighbors = (unsigned int)adjacencyListIndex.size();
        this->accessed = n2.accessed;
    }

    void swap(Node & other) noexcept {
        using std::swap;
        other.key = key;
        other.accessed = accessed;
        std::swap(adjacencyListIndex, other.adjacencyListIndex);
        other.totalNeighbors = totalNeighbors;
        other.rootNode = rootNode;
    }

    Node & operator=(const Node & nodeToCopy) {
        Node tmp(nodeToCopy);
        tmp.swap(*this);
        return *this;
    }


    Node(Node && other) noexcept  {
        key = other.key;
        totalNeighbors = other.totalNeighbors;
        rootNode = other.rootNode;
        accessed = other.accessed;
        adjacencyListIndex = std::move(other.adjacencyListIndex);
    }

    ~Node(){
        adjacencyListIndex.clear();
    }

    unsigned int getKey() const {return key;}
    void addNeighbor(Node * pNode);
    void removeNeighbor(Node * pNode);

    void printNeighbors();

    void setAccessed(bool value){this->accessed = value;}

    //void setPointerToTour(std::list < Node * > * pointer);//{ pointerToTour = pointer;} // this points to a map which should not invalidated
//    std::list < Node * > * getPointerToTour(){ return pointerToTour;}

    void setRootNodeOfTour(unsigned int index){rootNode = index;}
    unsigned int getRootNodeOfTour(){return rootNode;}

    bool getAccessed(){return this->accessed;}

    unsigned int getTotalNeighbors(){return adjacencyListIndex.size();}
//    Node * getPointerToNeighborByIndex(int index){return adjacencyList[index];}
    bool isNeighborPresent(unsigned int index);

    std::set<unsigned int > * getIteratorToIndices(){return &adjacencyListIndex;}

    bool validate();
};


#endif //IKETAMA_NODE_H
