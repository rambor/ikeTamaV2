//
// Created by xos81802 on 10/06/2021.
//
// Copyright 2021 Robert P. Rambo
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


#ifndef IKETAMA_DENSITYMAPPER_H
#define IKETAMA_DENSITYMAPPER_H

/*
 * DensityMapper: density map of an aligned set of bead models, refined against I(q) (and P(r)).
 *
 * The refinement itself is densrefine::Refiner (DensityRefiner.h): a continuous density on the bead lattice,
 * Metropolis single-point moves, I(q) from distance-shell sums updated incrementally. The earlier cross-entropy
 * search, the spherical-harmonic (SHE) intensity and the Debye-factor table have been removed (slower at these
 * model sizes; the shell sums are exact on the lattice).
 */

#include <string>
#include <vector>
#include <sastools/Datum.h>
#include <sastools/vector3.h>
#include <sastools/PofRData.h>
#include "DensityRefiner.h"

class DensityMapper {

    std::string inputFilename;                    // aligned set (EDGE RADIUS remark gives the bead lattice)
    std::vector<vector3> centered_coordinates;    // aligned set, centred

public:
    explicit DensityMapper(std::string filename);

    unsigned int getTotalCenteredCoordinates() const { return (unsigned int)centered_coordinates.size(); }

    /*
     * Density refinement against I(q) (and P(r)) starting from the occupancy of the aligned set (see DensityRefiner.h).
     *
     *   workingSet   I(q) to fit (q, I, sigma); averaged to options.pointsPerShannon points per Shannon channel
     *                (dmax from pofrData, else from the aligned set), alternate points held out as the free set
     *   pofrData     optional (nullptr): P(r) restraint, KL(P(r)_data || P(r)_model) on the Shannon bins
     *   options      massKDa is required; spacing <= 0 takes 2 x EDGE RADIUS of the aligned set
     *   prefix       writes <prefix>_map.xplor, _start_map.xplor, _density.pdb, _beads.pdb, _fit.dat
     *
     * Coordinates are in the centred frame of the aligned set.
     */
    densrefine::LambdaResult refineModelDensity(std::vector<Datum> & workingSet, PofRData * pofrData,
                                                densrefine::Options options, const std::string & prefix = "refined_density");

    static float readEdgeRadius(const std::string & pdbfile);
};

#endif //IKETAMA_DENSITYMAPPER_H
