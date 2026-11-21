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

#ifndef IKETAMA_SUBUNIT_H
#define IKETAMA_SUBUNIT_H

/*
 * SubUnit: a set of lattice beads (indices into the universe) with a rigid-body transform
 *
 *     transformed = R (x - pivot) + translateTo
 *
 * R is a 3x3 rotation (x1..z3, column-major as in the original code: new.x = x1*b.x + y1*b.y + z1*b.z).
 *
 * Changes from the 2020 version
 *  - copy/assignment are the compiler defaults: every member is copied, including the undo state (the hand-written
 *    versions left the backups empty, so undo after a copy restored zeros)
 *  - removeIndex, swapIndices, getIndexOfBeadtoRemove and getPointerToTransformedCoordinate throw std::logic_error
 *    when the index is not in the working set (they used to read/write past workingLimit)
 *  - capacity grows from 0 (an empty SubUnit could not take its first bead)
 *  - the constructor checks that indices has at least workingLimit entries
 *  - the random start uses the caller's generator (reproducible with --randomSeed) and is uniform over rotations
 *    (random unit quaternion; the old Euler draw in whole degrees over-sampled beta near +-90 deg)
 *  - rotateAt composes a small rotation with the current one (R <- dR R) instead of adding to Euler angles, which
 *    stalled near gimbal lock; inPlace = true keeps the subunit where it is when the pivot changes
 *  - backup/undo copy only the working set
 *  - translateToPosition(vector3) moves the pivot to an arbitrary point, e.g. a transformed bead of a neighbouring
 *    subunit (translateToTransformedPosition(index) uses the untransformed universe position of a bead, as before)
 *  - the declared-but-undefined constructor taking a translation vector is removed
 */

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <climits>
#include <cmath>
#include "sastools/vector3.h"
#include "sastools/Bead.h"
#include "PointSetModel.h"

class PointSetModel;

class SubUnit {

    std::vector<vector3> coordinates;                    // untransformed (universe) positions of the beads in use
    std::vector<vector3> transformedCoordinates;
    std::vector<vector3> backedUpTransformedCoordinates;

    // removing a coordinate swaps it to the end and decrements workingLimit; adding appends and increments
    std::vector<unsigned int> indicesInUse;

    unsigned int swapLocale = 0;

    unsigned int pivotIndex = 0, subUnitIndex = 0, translateTo = 0;
    unsigned int backedUpPivotIndex = 0, backedUpTranslateToIndex = 0;
    unsigned int capacity = 0;
    unsigned int workingLimit = 0;

    const Bead * pPivotPoint = nullptr, * pBackUpPivotPoint = nullptr;

    vector3 pivotVector, translateToVector, backedUpCenteringVector, backedUpTranslationVector;

    float x1 = 1, x2 = 0, x3 = 0, y1 = 0, y2 = 1, y3 = 0, z1 = 0, z2 = 0, z3 = 1;
    float x1_old = 1, x2_old = 0, x3_old = 0, y1_old = 0, y2_old = 1, y3_old = 0, z1_old = 0, z2_old = 0, z3_old = 1;

    static constexpr float convert = (float)(M_PI/180.0);

    void makeBackUp();
    void transformAll();                                  // recompute transformedCoordinates[0, workingLimit)
    void transformPoint(const vector3 *originalVec, vector3 *newVec) const;
    void setRotationMatrix(const float m[9]);             // row-major m -> x1..z3
    void getRotationMatrix(float m[9]) const;
    void ensureCapacity(unsigned int needed);
    unsigned int locate(unsigned int beadIndex, const char * caller) const;   // position in indicesInUse or throw

public:
    SubUnit() = default;

    SubUnit(unsigned int cofmIndex, unsigned int subUnitIndex, unsigned int translateTo);

    /**
     * @param indices    bead indices (at least workingLimit of them; the first workingLimit are used)
     * @param rng        generator for the random starting orientation (the run's generator, for reproducibility)
     */
    SubUnit(unsigned int subUnitIndex, unsigned int cofmIndex, unsigned int translateTo, unsigned int workingLimit,
            std::vector<unsigned int> &indices, PointSetModel *pModel, std::mt19937 &rng);

    SubUnit(const SubUnit &) = default;
    SubUnit & operator=(const SubUnit &) = default;
    SubUnit(SubUnit &&) noexcept = default;
    SubUnit & operator=(SubUnit &&) noexcept = default;
    ~SubUnit() = default;

    void addIndex(unsigned int indexToAdd, PointSetModel *pModel);
    void undoAdd();

    void removeIndex(unsigned int indexToRemove);
    void undoRemove();

    vector3 translateAndRotatePoint(const vector3 * originalVec) const;

    const vector3 *getPointerToTransformedCoordinates() const { return transformedCoordinates.data(); }
    const vector3 & getStartToTransformedCoordinates() const { return transformedCoordinates.front(); }

    unsigned int getIndexOfBeadtoRemove(unsigned int indexToRemove) const { return locate(indexToRemove, "getIndexOfBeadtoRemove"); }
    unsigned int getTranslateToIndex() const { return translateTo; }
    unsigned int getSubUnitIndex() const { return subUnitIndex; }
    unsigned int getWorkingLimit() const { return workingLimit; }
    unsigned int getPivotIndex() const { return pivotIndex; }
    const vector3 & getTranslateToVector() const { return translateToVector; }

    // pivot goes to the universe position of bead translatedTo (untransformed), as before
    void translateToTransformedPosition(unsigned int translatedTo, PointSetModel * pModel);
    // pivot goes to an arbitrary point (e.g. a transformed bead of a neighbouring subunit); translateTo index = UINT_MAX
    void translateToPosition(const vector3 & position);
    void undoTranslateTo(PointSetModel *pModel);

    /**
     * Rotate by small angles (degrees, about the lab x, y, z axes: dR = Rz Ry Rx) about bead newPivotPoint.
     * inPlace = false: the new pivot is placed at translateToVector (old behaviour: a pivot change also moves the subunit)
     * inPlace = true : the new pivot stays where it currently is, so the move is a pure rotation about that bead
     */
    void rotateAt(unsigned int newPivotPoint, float angleInDegreesX, float angleInDegreesY, float angleInDegreesZ,
                  PointSetModel * pModel, bool inPlace = false);
    void undoRotateTo();

    void swapIndices(unsigned int indexOfCurrent, unsigned int indexOfNew, PointSetModel * pModel);

    inline const vector3 & getPointerToTransformedCoordinate(unsigned int indexToFind) const {
        return transformedCoordinates[locate(indexToFind, "getPointerToTransformedCoordinate")];
    }
};

#endif //IKETAMA_SUBUNIT_H
