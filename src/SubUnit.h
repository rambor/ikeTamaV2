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

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "sastools/vector3.h"
#include "sastools/Bead.h"
#include "PointSetModel.h"
#include <random>

class PointSetModel;

class SubUnit {

    std::vector<vector3> backedUpCoordinates;
    std::vector<vector3> backedUpTransformedCoordinates;
    std::vector<vector3> coordinates;
    std::vector<vector3> transformedCoordinates;

    // removing a coordinate means swapping to the end and decrementing totalIndices
    // add means adding to totalIndices and incrementing

    std::vector<unsigned int> indicesInUse;
    std::vector<unsigned int> backedUpindicesInUse;

    unsigned int swapLocale;

    unsigned int pivotIndex, subUnitIndex, translateTo;
    unsigned int backedUpPivotIndex, backedUpTranslateToIndex;
    unsigned int capacity;

    Bead * pPivotPoint, * pBackUpPivotPoint;

    vector3 pivotVector, translateToVector, backedUpCenteringVector, backedUpTranslationVector;

    unsigned int workingLimit;
    float alpha_x, beta_y, gamma_z;
    float alpha_x_old, beta_y_old, gamma_z_old;
    float x1, x2, x3, y1, y2, y3, z1, z2, z3;
    float x1_old, x2_old, x3_old, y1_old, y2_old, y3_old, z1_old, z2_old, z3_old;

    float convert = (float)(M_PI/180.0f);

    void makeBackUp();

    void rotate(float alpha, float beta, float gamma);

    void transformPoint(const vector3 *originalVec, vector3 *newVewc);

public:
    SubUnit();

    SubUnit(unsigned int cofmIndex, unsigned int subUnitIndex, unsigned int translateTo);

    SubUnit(unsigned int subUnitIndex, unsigned int cofmIndex, unsigned translateTo, unsigned int workingLimit, std::vector<unsigned int> &indices, PointSetModel *pModel);

    SubUnit(unsigned int subUnitIndex, unsigned int cofmIndex, unsigned translateTo, vector3 * translationVec, unsigned int workingLimit, std::vector<unsigned int> &indices, PointSetModel *pModel);

    /**
     * Copy Constructor
     * @param p2
     */
    SubUnit(const SubUnit &p2) : pivotIndex(p2.pivotIndex), subUnitIndex(p2.subUnitIndex) {

        this->alpha_x = p2.alpha_x;
        this->beta_y = p2.beta_y;
        this->gamma_z = p2.gamma_z;

        this->capacity = p2.capacity;
        this->workingLimit = p2.workingLimit;

        this->x1 = p2.x1;
        this->x2 = p2.x2;
        this->x3 = p2.x3;
        this->y1 = p2.y1;
        this->y2 = p2.y2;
        this->y3 = p2.y3;
        this->z1 = p2.z1;
        this->z2 = p2.z2;
        this->z3 = p2.z3;

        this->pivotVector = vector3(p2.pivotVector.x, p2.pivotVector.y, p2.pivotVector.z);
        this->translateToVector = vector3(p2.translateToVector.x, p2.translateToVector.y, p2.translateToVector.z);

        this->pivotIndex = p2.pivotIndex;
        this->translateTo = p2.translateTo;

        this->pPivotPoint = *(&p2.pPivotPoint);
        this->pBackUpPivotPoint = *(&p2.pBackUpPivotPoint);

        this->coordinates.resize(this->capacity);
        this->transformedCoordinates.resize(this->capacity);
        this->indicesInUse.resize(this->capacity);
        this->backedUpCoordinates.resize(this->capacity);
        this->backedUpindicesInUse.resize(this->capacity);
        this->backedUpTransformedCoordinates.resize(this->capacity);

        std::copy(p2.coordinates.begin(),
                  p2.coordinates.begin() + this->workingLimit,
                  this->coordinates.begin());

        std::copy(p2.transformedCoordinates.begin(),
                  p2.transformedCoordinates.begin() + this->workingLimit,
                  this->transformedCoordinates.begin());

        std::copy(p2.indicesInUse.begin(),
                  p2.indicesInUse.begin() + this->workingLimit,
                  this->indicesInUse.begin());
    }


    // copy assignment
    SubUnit& operator=(const SubUnit& p2) {

        if (this == &p2)
            return *this;

        this->alpha_x = p2.alpha_x;
        this->beta_y = p2.beta_y;
        this->gamma_z = p2.gamma_z;

        this->capacity = p2.capacity;
        this->workingLimit = p2.workingLimit;

        this->x1 = p2.x1;
        this->x2 = p2.x2;
        this->x3 = p2.x3;
        this->y1 = p2.y1;
        this->y2 = p2.y2;
        this->y3 = p2.y3;
        this->z1 = p2.z1;
        this->z2 = p2.z2;
        this->z3 = p2.z3;

        this->pivotVector = vector3(p2.pivotVector.x, p2.pivotVector.y, p2.pivotVector.z);
        this->translateToVector = vector3(p2.translateToVector.x, p2.translateToVector.y, p2.translateToVector.z);

        this->pivotIndex = p2.pivotIndex;
        this->translateTo = p2.translateTo;

        this->pPivotPoint = *(&p2.pPivotPoint);
        this->pBackUpPivotPoint = *(&p2.pBackUpPivotPoint);

        this->coordinates.resize(this->capacity);
        this->transformedCoordinates.resize(this->capacity);
        this->indicesInUse.resize(this->capacity);
        this->backedUpCoordinates.resize(this->capacity);
        this->backedUpindicesInUse.resize(this->capacity);
        this->backedUpTransformedCoordinates.resize(this->capacity);

        std::copy(p2.coordinates.begin(),
                  p2.coordinates.begin() + this->workingLimit,
                  this->coordinates.begin());

        std::copy(p2.transformedCoordinates.begin(),
                  p2.transformedCoordinates.begin() + this->workingLimit,
                  this->transformedCoordinates.begin());

        std::copy(p2.indicesInUse.begin(),
                  p2.indicesInUse.begin() + this->workingLimit,
                  this->indicesInUse.begin());

        return *this;
    }


    ~SubUnit() = default;


    void addIndex(unsigned int indexToAdd, PointSetModel *pModel);

    void undoAdd();

    void removeIndex(unsigned int indexToRemove);

    void undoRemove();

    vector3 translateAndRotatePoint(const vector3 * const originalVec);

    const vector3 *getPointerToTransformedCoordinates() const { return &transformedCoordinates.front(); }

    unsigned int getIndexOfBeadtoRemove(unsigned int indexToRemove);
    unsigned int getTranslateToIndex(){ return translateTo;}
    unsigned int getSubUnitIndex() { return subUnitIndex; }

    void translateToTransformedPosition(unsigned int translatedTo, PointSetModel * pModel);

    void undoTranslateTo(PointSetModel *pModel);

    void rotateAt(unsigned int newPivotPoint, float angleInDegreesX, float angleInDegreesY, float angleInDegreesZ, PointSetModel * pModel);

    void undoRotateTo();
    void swapIndices(unsigned int indexOfCurrent, unsigned int indexOfNew, PointSetModel * pModel);

    const vector3 & getStartToTransformedCoordinates() const {return transformedCoordinates.front();}

    inline vector3 & getPointerToTransformedCoordinate(int indexToFind){
        auto locale = std::find(indicesInUse.begin(), indicesInUse.begin() + workingLimit, indexToFind);
        int index = std::distance(indicesInUse.begin(), locale);

//        if (index >= workingLimit){
//            auto actual = std::find(indicesInUse.begin(), indicesInUse.end(), indexToFind);
//            std::cout << " CAN NOT FIND  " << indexToFind << " => "<< index << " " << workingLimit << std::endl;
//            std::cout << "    DISTANCE   " << std::distance(indicesInUse.begin(), actual) << " => " << " " << workingLimit << " | " << indicesInUse.size() << std::endl;
//            exit(0);
//        }
        return transformedCoordinates[index];
    };
};


#endif //IKETAMA_SUBUNIT_H
