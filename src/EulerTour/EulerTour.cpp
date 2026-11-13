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


#include "EulerTour.h"
#include "../PointSetModel.h"
#include <algorithm>

// Constructor
EulerTour::EulerTour() = default;

EulerTour::EulerTour(std::vector<unsigned int> & indices, const unsigned int workingLimit, PointSetModel *pModel){
    this->createInitialTour(workingLimit, pModel, indices);
}

EulerTour::EulerTour(std::vector<unsigned int>::const_iterator beginIt, const unsigned int workingLimit, PointSetModel *pModel){
    attach(pModel);
    labelAll(std::vector<unsigned int>(beginIt, beginIt + workingLimit));
}

EulerTour::EulerTour(std::set<unsigned int> & indices, PointSetModel *pModel){
    attach(pModel);
    labelAll(std::vector<unsigned int>(indices.begin(), indices.end()));
}

void EulerTour::createInitialTour(unsigned int workingLimit, PointSetModel *pModel, std::vector<unsigned int> & indices) {
    attach(pModel);
    labelAll(std::vector<unsigned int>(indices.begin(), indices.begin() + workingLimit));
}

unsigned int EulerTour::newTour(std::vector<unsigned int>::iterator beginIt, unsigned int workingLimit, PointSetModel *pModel) {
    attach(pModel);
    labelAll(std::vector<unsigned int>(beginIt, beginIt + workingLimit));
    return totalComponents;
}


void EulerTour::attach(PointSetModel * model){
    pModel = model;
    neighborhood = model->getDirectPointerToNeighborhood();
    sizeOfNeighborhood = model->getSizeOfNeighborhood();
    neighborLimit = model->getNeighborLimit();
    universe = model->getTotalNumberOfBeadsInUniverse();
    label.assign(universe, -1);
    stamp.assign(universe, 0);
    owner.assign(universe, 0);
    compSize.clear(); freeLabels.clear();
    totalComponents = 0; totalNodes = 0; epoch = 0;
}


int EulerTour::newLabel(){
    if (!freeLabels.empty()){
        int l = freeLabels.back(); freeLabels.pop_back(); compSize[l] = 0; return l;
    }
    compSize.push_back(0);
    return (int)compSize.size() - 1;
}

unsigned int EulerTour::nextEpoch(){
    if (++epoch == 0){ std::fill(stamp.begin(), stamp.end(), 0); epoch = 1; }
    return epoch;
}


/**
 * labels every point of the set by breadth-first search
 */
void EulerTour::labelAll(const std::vector<unsigned int> & points){
    for (unsigned int p : points) label[p] = -2;   // present, not yet labelled
    totalNodes = (unsigned int)points.size();
    std::vector<unsigned int> queue;
    queue.reserve(points.size());
    for (unsigned int p : points){
        if (label[p] != -2) continue;
        const int l = newLabel();
        queue.clear(); queue.push_back(p); label[p] = l;
        for (size_t h = 0; h < queue.size(); h++){
            const unsigned int * nb = neighborhood + (size_t)sizeOfNeighborhood*queue[h];
            for (unsigned int k = 0; k < sizeOfNeighborhood; k++){
                const unsigned int w = nb[k];
                if (w == neighborLimit) break;
                if (w < neighborLimit && label[w] == -2){ label[w] = l; queue.push_back(w); }
            }
        }
        compSize[l] = (unsigned int)queue.size();
        totalComponents++;
    }
}


void EulerTour::relabelComponent(unsigned int start, int from, int to){
    std::vector<unsigned int> queue;
    queue.reserve(compSize[from]);
    queue.push_back(start); label[start] = to;
    for (size_t h = 0; h < queue.size(); h++){
        const unsigned int * nb = neighborhood + (size_t)sizeOfNeighborhood*queue[h];
        for (unsigned int k = 0; k < sizeOfNeighborhood; k++){
            const unsigned int w = nb[k];
            if (w == neighborLimit) break;
            if (w < neighborLimit && label[w] == from){ label[w] = to; queue.push_back(w); }
        }
    }
}


/**
 * add a node to existing set of nodes, merge the components it bridges
 * returns number of components
 */
unsigned int EulerTour::addNode(unsigned int newNode, PointSetModel * model) {

    if (pModel == nullptr) attach(model);
    if (label[newNode] >= 0) return totalComponents;   // already in the set

    // distinct components among the occupied neighbours (at most 12)
    int labels[64]; unsigned int seeds[64]; unsigned int distinct = 0;
    const unsigned int * nb = neighborhood + (size_t)sizeOfNeighborhood*newNode;
    for (unsigned int k = 0; k < sizeOfNeighborhood; k++){
        const unsigned int w = nb[k];
        if (w == neighborLimit) break;
        if (w < neighborLimit && label[w] >= 0){
            bool seen = false;
            for (unsigned int d = 0; d < distinct; d++) if (labels[d] == label[w]) { seen = true; break; }
            if (!seen){ labels[distinct] = label[w]; seeds[distinct] = w; distinct++; }
        }
    }

    totalNodes++;
    if (distinct == 0){
        const int l = newLabel();
        label[newNode] = l; compSize[l] = 1;
        totalComponents++;
        return totalComponents;
    }

    // largest neighbouring component absorbs the others
    unsigned int best = 0;
    for (unsigned int d = 1; d < distinct; d++) if (compSize[labels[d]] > compSize[labels[best]]) best = d;
    const int target = labels[best];
    label[newNode] = target; compSize[target]++;

    for (unsigned int d = 0; d < distinct; d++){
        if (d == best) continue;
        relabelComponent(seeds[d], labels[d], target);
        compSize[target] += compSize[labels[d]];
        compSize[labels[d]] = 0;
        freeLabels.push_back(labels[d]);
        totalComponents--;
    }
    return totalComponents;
}


/**
 * @param indexOfNode
 * @return total number of components after removal
 */
unsigned int EulerTour::removeNode(unsigned int indexOfNode){

    const int L = label[indexOfNode];
    if (L < 0) return totalComponents;   // not in the set
    label[indexOfNode] = -1;
    compSize[L]--;
    totalNodes--;

    // occupied neighbours
    unsigned int nbrs[64]; unsigned int k = 0;
    const unsigned int * nb = neighborhood + (size_t)sizeOfNeighborhood*indexOfNode;
    for (unsigned int j = 0; j < sizeOfNeighborhood; j++){
        const unsigned int w = nb[j];
        if (w == neighborLimit) break;
        if (w < neighborLimit && label[w] >= 0) nbrs[k++] = w;
    }

    if (k == 0){                                   // isolated point: its component disappears
        freeLabels.push_back(L);
        totalComponents--;
        return totalComponents;
    }
    if (k == 1) return totalComponents;            // a leaf: nothing can split

    // 1. local test: are the neighbours connected to each other directly (within the neighbourhood)?
    unsigned int parent[64];
    for (unsigned int a = 0; a < k; a++) parent[a] = a;
    auto findRoot = [&](unsigned int a){ while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; };
    unsigned int groups = k;
    for (unsigned int a = 0; a < k; a++){
        const unsigned int * na = neighborhood + (size_t)sizeOfNeighborhood*nbrs[a];
        for (unsigned int j = 0; j < sizeOfNeighborhood; j++){
            const unsigned int w = na[j];
            if (w == neighborLimit) break;
            for (unsigned int b = a + 1; b < k; b++){
                if (nbrs[b] == w){
                    const unsigned int ra = findRoot(a), rb = findRoot(b);
                    if (ra != rb){ parent[rb] = ra; groups--; }
                }
            }
        }
    }
    if (groups == 1) return totalComponents;

    // 2. one breadth-first search per local group, expanded in turn; searches that meet are merged.
    //    Stop when all searches belong to one class (still connected) or when at most one class is unfinished.
    unsigned int groupOf[64], G = 0;
    {
        unsigned int rootToGroup[64];
        for (unsigned int a = 0; a < k; a++) rootToGroup[a] = 64;
        for (unsigned int a = 0; a < k; a++){
            const unsigned int r = findRoot(a);
            if (rootToGroup[r] == 64) rootToGroup[r] = G++;
            groupOf[a] = rootToGroup[r];
        }
    }

    const unsigned int e = nextEpoch();
    std::vector<std::vector<unsigned int>> queue(G);
    std::vector<size_t> head(G, 0);
    std::vector<unsigned int> cls(G);
    for (unsigned int g = 0; g < G; g++) cls[g] = g;
    auto classOf = [&](unsigned int g){ while (cls[g] != g) g = cls[g] = cls[cls[g]]; return g; };

    for (unsigned int a = 0; a < k; a++){
        stamp[nbrs[a]] = e; owner[nbrs[a]] = (unsigned short)groupOf[a];
        queue[groupOf[a]].push_back(nbrs[a]);
    }

    unsigned int classes = G;
    std::vector<char> classActive(G);
    while (true){
        // which classes still have work
        std::fill(classActive.begin(), classActive.end(), 0);
        unsigned int unfinished = 0;
        for (unsigned int g = 0; g < G; g++){
            if (head[g] < queue[g].size()){
                const unsigned int c = classOf(g);
                if (!classActive[c]){ classActive[c] = 1; unfinished++; }
            }
        }
        if (unfinished <= 1) break;

        for (unsigned int g = 0; g < G; g++){        // one expansion per search
            if (head[g] >= queue[g].size()) continue;
            const unsigned int u = queue[g][head[g]++];
            const unsigned int * nu = neighborhood + (size_t)sizeOfNeighborhood*u;
            for (unsigned int j = 0; j < sizeOfNeighborhood; j++){
                const unsigned int w = nu[j];
                if (w == neighborLimit) break;
                if (w >= neighborLimit || label[w] < 0) continue;
                if (stamp[w] != e){
                    stamp[w] = e; owner[w] = (unsigned short)g; queue[g].push_back(w);
                } else {
                    const unsigned int c1 = classOf(g), c2 = classOf(owner[w]);
                    if (c1 != c2){
                        cls[c2] = c1;
                        if (--classes == 1) return totalComponents;   // all pieces meet: still connected
                    }
                }
            }
        }
    }

    // split: every class except the unfinished one (or, if all finished, the largest) is a separate component
    std::vector<unsigned int> classSize(G, 0);
    std::vector<char> unfinishedClass(G, 0);
    for (unsigned int g = 0; g < G; g++){
        const unsigned int c = classOf(g);
        classSize[c] += (unsigned int)queue[g].size();
        if (head[g] < queue[g].size()) unfinishedClass[c] = 1;
    }
    int keep = -1;
    for (unsigned int c = 0; c < G; c++) if (classOf(c) == c && unfinishedClass[c]) keep = (int)c;
    if (keep < 0){
        for (unsigned int c = 0; c < G; c++)
            if (classOf(c) == c && (keep < 0 || classSize[c] > classSize[keep])) keep = (int)c;
    }

    std::vector<int> newLabelOfClass(G, -1);
    for (unsigned int c = 0; c < G; c++){
        if (classOf(c) != c || (int)c == keep) continue;
        const int l = newLabel();
        newLabelOfClass[c] = l;
        compSize[l] = classSize[c];
        compSize[L] -= classSize[c];
        totalComponents++;
    }
    for (unsigned int g = 0; g < G; g++){
        const int l = newLabelOfClass[classOf(g)];
        if (l < 0) continue;
        for (unsigned int u : queue[g]) label[u] = l;
    }
    return totalComponents;
}


unsigned int EulerTour::countComponentsBruteForce() const {
    std::vector<char> seen(universe, 0);
    std::vector<unsigned int> queue;
    unsigned int count = 0;
    for (unsigned int p = 0; p < universe; p++){
        if (label[p] < 0 || seen[p]) continue;
        count++;
        queue.clear(); queue.push_back(p); seen[p] = 1;
        for (size_t h = 0; h < queue.size(); h++){
            const unsigned int * nb = neighborhood + (size_t)sizeOfNeighborhood*queue[h];
            for (unsigned int k = 0; k < sizeOfNeighborhood; k++){
                const unsigned int w = nb[k];
                if (w == neighborLimit) break;
                if (w < neighborLimit && label[w] >= 0 && !seen[w]){ seen[w] = 1; queue.push_back(w); }
            }
        }
    }
    return count;
}


bool EulerTour::validateList(std::string text){
    const unsigned int count = countComponentsBruteForce();
    if (count != totalComponents){
        std::cout << " @ " << text << " COMPONENTS " << totalComponents << " != RECOUNT " << count << std::endl;
        return true;
    }
    return false;
}

bool EulerTour::validateNodes(std::string str){
    return !validateList(str);
}


void EulerTour::buildCompatibilityViews(){
    tours.clear();
    nodes.clear();
    for (unsigned int p = 0; p < universe; p++){
        if (label[p] < 0) continue;
        Node & node = nodes.emplace(p, Node(p)).first->second;
        node.setRootNodeOfTour((unsigned int)label[p]);
    }
    for (auto & entry : nodes){
        const unsigned int * nb = neighborhood + (size_t)sizeOfNeighborhood*entry.first;
        for (unsigned int k = 0; k < sizeOfNeighborhood; k++){
            const unsigned int w = nb[k];
            if (w == neighborLimit) break;
            auto it = nodes.find(w);
            if (w < neighborLimit && it != nodes.end()) entry.second.addNeighbor(&it->second);
        }
    }
    // one list per component, keyed by its smallest lattice index
    std::map<int, unsigned int> rootOfLabel;
    for (auto & entry : nodes){
        const int l = label[entry.first];
        auto r = rootOfLabel.emplace(l, entry.first).first;
        tours[r->second].push_back(&entry.second);
        entry.second.setRootNodeOfTour(r->second);
    }
}
