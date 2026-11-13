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
// Include after Anneal.h (needs qhull and PointSetModel).
#ifndef IKETAMA_HULLVOLUMECACHE_H
#define IKETAMA_HULLVOLUMECACHE_H

#include "Anneal.h"

/*
 * Convex hull volume of the selected points, recomputed only when the hull can change.
 * When qhull is called it is called exactly as calculateCVXHULLVolume does (same points, same order, same flags).
 */
class HullVolumeCache {
    struct State {
        float volume = 0.0f;
        std::vector<double> planes;          // nx, ny, nz, offset per facet (outward normals)
        std::vector<unsigned int> vertices;  // lattice indices of hull vertices
    };
    State committed, pending;
    bool hasPending = false;
    std::vector<uint8_t> isVertex;
    std::vector<coordT> points;
    PointSetModel * pModel;
    char * flags;
    static constexpr double eps = 1.0e-4;    // Angstrom; points this close to a facet count as inside

    void compute(const std::vector<unsigned int> & v, unsigned int wl, State & out){
        points.resize(3*(size_t)wl);
        for (unsigned int i=0; i<wl; i++){
            Bead * pBead = pModel->getBead(v[i]);
            points[3*i] = pBead->getX(); points[3*i+1] = pBead->getY(); points[3*i+2] = pBead->getZ();
        }
        qh_new_qhull(3, (int)wl, points.data(), 0, flags, nullptr, nullptr);
        out.volume = (float)(qh totvol);
        out.planes.clear(); out.vertices.clear();
        facetT * facet;
        FORALLfacets {
            out.planes.push_back(facet->normal[0]); out.planes.push_back(facet->normal[1]);
            out.planes.push_back(facet->normal[2]); out.planes.push_back(facet->offset);
        }
        vertexT * vertex;
        FORALLvertices {
            out.vertices.push_back(v[qh_pointid(vertex->point)]);
        }
        qh_freeqhull(true);
    }

    bool inside(unsigned int x) const {
        Bead * pBead = pModel->getBead(x);
        const double px = pBead->getX(), py = pBead->getY(), pz = pBead->getZ();
        const double * p = committed.planes.data();
        for (size_t f = 0; f < committed.planes.size(); f += 4){
            if (p[f]*px + p[f+1]*py + p[f+2]*pz + p[f+3] > eps) return false;
        }
        return true;
    }

    float recompute(const std::vector<unsigned int> & v, unsigned int wl){
        compute(v, wl, pending);
        hasPending = true;
        return pending.volume;
    }

public:
    unsigned int recomputed = 0, reused = 0;

    HullVolumeCache(PointSetModel * model, char * qhullFlags, const std::vector<unsigned int> & v, unsigned int wl)
            : isVertex(model->getTotalNumberOfBeadsInUniverse(), 0), pModel(model), flags(qhullFlags) {
        compute(v, wl, committed);
        for (unsigned int i : committed.vertices) isVertex[i] = 1;
    }

    float currentVolume() const { return committed.volume; }

    // v/wl already describe the trial set
    float afterAdd(unsigned int x, const std::vector<unsigned int> & v, unsigned int wl){
        if (inside(x)) { hasPending = false; reused++; return committed.volume; }
        recomputed++; return recompute(v, wl);
    }
    float afterRemove(unsigned int x, const std::vector<unsigned int> & v, unsigned int wl){
        if (!isVertex[x]) { hasPending = false; reused++; return committed.volume; }
        recomputed++; return recompute(v, wl);
    }
    float afterSwap(unsigned int out, unsigned int in, const std::vector<unsigned int> & v, unsigned int wl){
        if (!isVertex[out] && inside(in)) { hasPending = false; reused++; return committed.volume; }
        recomputed++; return recompute(v, wl);
    }
    void accept(){
        if (hasPending){
            for (unsigned int i : committed.vertices) isVertex[i] = 0;
            std::swap(committed, pending);
            for (unsigned int i : committed.vertices) isVertex[i] = 1;
        }
        hasPending = false;
    }
    void reject(){ hasPending = false; }
};


#endif //IKETAMA_HULLVOLUMECACHE_H
