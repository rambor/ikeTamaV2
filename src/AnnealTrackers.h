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
// Incremental bookkeeping used by the annealing loops. Each class replaces a full recomputation per move with
// an update of the moved bead only, and is exact: it returns what the full recomputation would.
//
//   ContactTracker   occupancy + contact number of every lattice point, contact energy
//                    replaces std::set beads_in_use_tree, numberOfContactsFromSet, calculateNumberOfPairwiseContacts
//   SortedSelection  the sorted bead_indices vector with O(1) location of free points and exact undo
//                    replaces std::find over the free points, std::sort after every move and the backUpState copies
//   SwapSelection    same for the symmetry code, where bead_indices is not kept sorted
//   ClashTracker     number of clashes between the asymmetric unit and its symmetry mates
//                    replaces getViolations / getViolationsFromSet (O(S*N^2) per move)
//
#ifndef IKETAMA_ANNEALTRACKERS_H
#define IKETAMA_ANNEALTRACKERS_H

#include <array>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <unordered_map>

/*
 * Contact potential per bead as a function of its number of occupied neighbours.
 * calculateNumberOfPairwiseContacts (full recount) and ContactTracker (incremental) both use this,
 * so edit the potential here only.
 */
inline double contactCost(unsigned int val){
    double diff = 0;
    if (val < 4) {
        diff = 4.5 - val;
    } else if (val > 9.5) {
        diff = val - 9.0;
    }
    if (val > 3) {               // overrides the line above for val > 3 (as in the original function)
        diff = val - 3.0;
    }
    return diff*diff;
}


class ContactTracker {
    const unsigned int * nbr = nullptr;  // neighbour table nbr[K*i + k], list ends at 'limit'
    unsigned int K = 12, limit = 0, universe = 0;
    std::vector<uint8_t> occ;            // 1 if lattice point is in the model
    std::vector<uint8_t> cnt;            // number of occupied neighbours, kept for EVERY lattice point
    std::array<double, 13> cost{};
    std::array<unsigned int, 13> hist{}; // number of occupied points with c occupied neighbours
    unsigned int n = 0;

public:
    ContactTracker() = default;
    ContactTracker(const unsigned int * neighbours, unsigned int sizeOfNeighbourhood, unsigned int neighbourLimit, unsigned int totalPoints)
            : nbr(neighbours), K(sizeOfNeighbourhood), limit(neighbourLimit), universe(totalPoints),
              occ(totalPoints, 0), cnt(totalPoints, 0) {
        for (unsigned int c = 0; c < 13; c++) cost[c] = contactCost(c);
    }

    void clear(){ std::fill(occ.begin(), occ.end(), 0); std::fill(cnt.begin(), cnt.end(), 0); hist.fill(0); n = 0; }

    template<class It> void reset(It begin, unsigned int workingLimit){
        clear();
        for (unsigned int b = 0; b < workingLimit; b++) add(*(begin + b));
    }

    inline bool occupied(unsigned int i) const { return occ[i] != 0; }
    inline unsigned int contacts(unsigned int i) const { return cnt[i]; }   // occupied neighbours of any point
    inline unsigned int size() const { return n; }
    inline const std::array<unsigned int, 13> & histogram() const { return hist; }

    // mean contact cost over occupied beads == calculateNumberOfPairwiseContacts (exactly, for the default potential)
    inline double energy() const {
        double sum = 0.0;
        for (unsigned int c = 0; c < 13; c++) sum += hist[c]*cost[c];
        return n ? sum/(double)n : 0.0;
    }

    // free, valid neighbours of i (== number counted by getProbabilityOfSelectedIndex)
    inline unsigned int freeNeighbours(unsigned int i) const {
        unsigned int valid = 0;
        const unsigned int * p = nbr + (size_t)K*i;
        for (unsigned int k = 0; k < K; k++){
            if (p[k] == limit) break;
            if (p[k] < limit) valid++;
        }
        return valid - cnt[i];
    }

    inline void add(unsigned int i){
        assert(!occ[i]);
        const unsigned int * p = nbr + (size_t)K*i;
        for (unsigned int k = 0; k < K; k++){
            const unsigned int j = p[k];
            if (j == limit) break;
            if (j > limit) continue;
            if (occ[j]) { hist[cnt[j]]--; hist[cnt[j] + 1]++; }
            cnt[j]++;
        }
        occ[i] = 1; n++;
        hist[cnt[i]]++;
    }

    inline void remove(unsigned int i){
        assert(occ[i]);
        hist[cnt[i]]--;
        occ[i] = 0; n--;
        const unsigned int * p = nbr + (size_t)K*i;
        for (unsigned int k = 0; k < K; k++){
            const unsigned int j = p[k];
            if (j == limit) break;
            if (j > limit) continue;
            cnt[j]--;
            if (occ[j]) { hist[cnt[j] + 1]--; hist[cnt[j]]++; }
        }
    }
};


/*
 * bead_indices kept as in the original code: [0, workingLimit) sorted, the rest in arbitrary order.
 * Every operation leaves the vector in exactly the state the original code (iter_swap + std::sort, restore from
 * backUpState) produces, so all later random draws pick the same beads. Cost O(workingLimit) memmove instead of
 * O(N_universe) find + copy and O(n log n) sort.
 */
class SortedSelection {
    std::vector<unsigned int> & v;
    unsigned int & wl;
    std::vector<unsigned int> tailPos;  // position in v, valid only for points at or beyond wl
    unsigned int lastP = 0, lastY = 0, lastQ = 0;

    inline void insertSorted(unsigned int x){             // into [0,wl), grows wl
        const unsigned int k = find(x);
        std::memmove(&v[k + 1], &v[k], (size_t)(wl - k)*sizeof(unsigned int));
        v[k] = x; wl++;
    }
    inline void eraseSorted(unsigned int k){              // shrinks wl
        std::memmove(&v[k], &v[k + 1], (size_t)(wl - k - 1)*sizeof(unsigned int));
        wl--;
    }

public:
    SortedSelection(std::vector<unsigned int> & vec, unsigned int & workingLimit, unsigned int universe)
            : v(vec), wl(workingLimit), tailPos(universe, 0) {
        for (unsigned int i = wl; i < (unsigned int)v.size(); i++) tailPos[v[i]] = i;
    }

    inline unsigned int find(unsigned int x) const {      // index of x in the sorted prefix
        return (unsigned int)(std::lower_bound(v.begin(), v.begin() + wl, x) - v.begin());
    }

    // == iter_swap(v[wl], find(v.begin()+wl, v.end(), x)); wl++; sort(v, v+wl)
    inline void add(unsigned int x){
        const unsigned int p = tailPos[x], y = v[wl];
        v[p] = y; tailPos[y] = p;
        lastP = p; lastY = y;
        insertSorted(x);
    }
    // exact inverse of the preceding add(x)  (== restoreAddingFromBackUp)
    inline void undoAdd(unsigned int x){
        eraseSorted(find(x));
        v[wl] = lastY; v[lastP] = x;
        tailPos[lastY] = wl; tailPos[x] = lastP;
    }

    // == wl--; iter_swap(v[k], v[wl]); sort(v, v+wl)      (removeLatticePositionToModel / ...ByIndex)
    inline void removeAt(unsigned int k){
        const unsigned int x = v[k];
        eraseSorted(k);
        v[wl] = x; tailPos[x] = wl;
    }
    // == wl++; sort(v, v+wl)                              (restoreRemovingLatticePointFromBackUp)
    inline void undoRemove(){ insertSorted(v[wl]); }

    // x = v[k] leaves, y (free) enters: == iter_swap(v[k], find(tail, y)); sort(v, v+wl)
    inline void swapAt(unsigned int k, unsigned int y){
        const unsigned int x = v[k], q = tailPos[y];
        v[q] = x; tailPos[x] = q; lastQ = q;
        eraseSorted(k);
        insertSorted(y);
    }
    // exact inverse of the preceding swapAt (== copy from backUpState)
    inline void undoSwap(unsigned int x, unsigned int y){
        eraseSorted(find(y));
        insertSorted(x);
        v[lastQ] = y; tailPos[y] = lastQ;
    }
};


/*
 * bead_indices for the symmetry code, where the prefix is not sorted (only swaps).
 * Keeps the position of every lattice point, so no linear searches are needed.
 */
class SwapSelection {
    std::vector<unsigned int> & v;
    unsigned int & wl;
    std::vector<unsigned int> pos;      // position of every lattice point in v
    unsigned int lastA = 0, lastB = 0;

    inline void swapPositions(unsigned int a, unsigned int b){
        std::swap(v[a], v[b]);
        pos[v[a]] = a; pos[v[b]] = b;
    }

public:
    SwapSelection(std::vector<unsigned int> & vec, unsigned int & workingLimit, unsigned int universe)
            : v(vec), wl(workingLimit), pos(universe, 0) {
        for (unsigned int i = 0; i < (unsigned int)v.size(); i++) pos[v[i]] = i;
    }
    inline unsigned int position(unsigned int x) const { return pos[x]; }

    // == iter_swap(v[wl], find(tail, x)); wl++
    inline void add(unsigned int x){ lastA = wl; lastB = pos[x]; swapPositions(wl, pos[x]); wl++; }
    inline void undoAdd(){ wl--; swapPositions(lastA, lastB); }

    // == wl--; iter_swap(find(prefix, x), v[wl])         (removeLatticePositionToModelSym)
    inline void remove(unsigned int x){ wl--; lastA = pos[x]; lastB = wl; swapPositions(pos[x], wl); }
    inline void undoRemove(){ swapPositions(lastA, lastB); wl++; }

    // == iter_swap(v[k], find(tail, y))
    inline void swapAt(unsigned int k, unsigned int y){ lastA = k; lastB = pos[y]; swapPositions(k, pos[y]); }
    inline void undoSwap(){ swapPositions(lastA, lastB); }
};


/*
 * Clashes between the asymmetric unit and its symmetry mates, as counted by getViolations:
 *     V = sum over beads i, mates s >= 1, beads j of [ |b_i - R_s b_j| < limit ]
 * The lattice and the symmetry operators are fixed, so the partners of every lattice point are found once
 * (with the same coordinates and comparison as getViolations) and each move updates V from the moved bead only.
 */
class ClashTracker {
    std::vector<unsigned int> fStart, fList;  // forward:  j in F(i)  <=> |b_i - R_s b_j| < limit for some s (with multiplicity)
    std::vector<unsigned int> rStart, rList;  // reverse:  i in R(j)  <=> j in F(i)
    std::vector<uint8_t> occ;
    unsigned int violations = 0;
    bool built = false;

    inline unsigned int delta(unsigned int a) const {   // pairs involving a, with a occupied
        unsigned int d = 0;
        for (unsigned int t = fStart[a]; t < fStart[a + 1]; t++) d += occ[fList[t]];
        for (unsigned int t = rStart[a]; t < rStart[a + 1]; t++) { const unsigned int i = rList[t]; if (i != a) d += occ[i]; }
        return d;
    }

public:
    bool isBuilt() const { return built; }
    unsigned int count() const { return violations; }

    /*
     * coords: lattice coordinates (N); mates: coordinates of R_s b_j for s = 1..S-1, laid out as mates[(s-1)*N + j];
     * both produced by PointSetModel exactly as getViolations does. VEC must provide x, y, z and (a - b).length().
     */
    template<class VEC>
    void build(const std::vector<VEC> & coords, const std::vector<VEC> & mates, unsigned int S, float limit, const std::vector<char> * useMate = nullptr){
        const unsigned int N = (unsigned int)coords.size();
        // spatial hash of the lattice points, cell size = limit
        const float cell = std::max(limit, 1e-3f), inv = 1.0f/cell;
        auto key = [&](int x, int y, int z){ return ((int64_t)(x + 1048576) << 42) | ((int64_t)(y + 1048576) << 21) | (int64_t)(z + 1048576); };
        std::unordered_map<int64_t, std::vector<unsigned int>> grid;
        grid.reserve(N);
        for (unsigned int i = 0; i < N; i++)
            grid[key((int)std::floor(coords[i].x*inv), (int)std::floor(coords[i].y*inv), (int)std::floor(coords[i].z*inv))].push_back(i);

        std::vector<std::vector<unsigned int>> F(N);
        for (unsigned int s = 1; s < S; s++){
            if (useMate && (s >= useMate->size() || !(*useMate)[s])) continue;   // only the selected mates (interfaces)
            for (unsigned int j = 0; j < N; j++){
                const VEC & m = mates[(size_t)(s - 1)*N + j];
                const int cx = (int)std::floor(m.x*inv), cy = (int)std::floor(m.y*inv), cz = (int)std::floor(m.z*inv);
                for (int dx = -1; dx <= 1; dx++) for (int dy = -1; dy <= 1; dy++) for (int dz = -1; dz <= 1; dz++){
                    auto it = grid.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == grid.end()) continue;
                    for (unsigned int i : it->second){
                        const float distance_to = (coords[i] - m).length();   // same expression as getViolations
                        if (distance_to < limit) F[i].push_back(j);
                    }
                }
            }
        }
        fStart.assign(N + 1, 0); rStart.assign(N + 1, 0);
        for (unsigned int i = 0; i < N; i++){ fStart[i + 1] = fStart[i] + (unsigned int)F[i].size(); for (unsigned int j : F[i]) rStart[j + 1]++; }
        for (unsigned int j = 0; j < N; j++) rStart[j + 1] += rStart[j];
        fList.resize(fStart[N]); rList.resize(rStart[N]);
        std::vector<unsigned int> fill(rStart.begin(), rStart.end() - 1);
        for (unsigned int i = 0; i < N; i++){
            std::copy(F[i].begin(), F[i].end(), fList.begin() + fStart[i]);
            for (unsigned int j : F[i]) rList[fill[j]++] = i;
        }
        occ.assign(N, 0); violations = 0; built = true;
    }

    template<class It> void reset(It begin, unsigned int workingLimit){
        std::fill(occ.begin(), occ.end(), 0); violations = 0;
        for (unsigned int b = 0; b < workingLimit; b++) add(*(begin + b));
    }
    // clashes of an arbitrary set (same value as getViolations), without changing the tracked state
    template<class It> unsigned int countFor(It begin, unsigned int workingLimit){
        std::vector<uint8_t> saved(occ);
        std::fill(occ.begin(), occ.end(), 0);
        It it = begin;
        for (unsigned int b = 0; b < workingLimit; b++, ++it) occ[*it] = 1;
        unsigned int total = 0;
        it = begin;
        for (unsigned int b = 0; b < workingLimit; b++, ++it){
            const unsigned int i = *it;
            for (unsigned int t = fStart[i]; t < fStart[i + 1]; t++) total += occ[fList[t]];
        }
        occ.swap(saved);
        return total;
    }
    // thread-safe count for a set given as lattice indices; 'mark' is caller-owned scratch (size N, all zero, left zero)
    unsigned int countWith(std::vector<uint8_t> & mark, const unsigned int * beads, unsigned int n) const {
        for (unsigned int b = 0; b < n; b++) mark[beads[b]] = 1;
        unsigned int total = 0;
        for (unsigned int b = 0; b < n; b++){
            const unsigned int i = beads[b];
            for (unsigned int t = fStart[i]; t < fStart[i + 1]; t++) total += mark[fList[t]];
        }
        for (unsigned int b = 0; b < n; b++) mark[beads[b]] = 0;
        return total;
    }
    inline void add(unsigned int a){ occ[a] = 1; violations += delta(a); }
    inline void remove(unsigned int a){ violations -= delta(a); occ[a] = 0; }
};

#endif //IKETAMA_ANNEALTRACKERS_H
