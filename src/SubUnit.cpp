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

#include "SubUnit.h"

SubUnit::SubUnit(unsigned int pivotIndex, unsigned int subUnitIndex, unsigned int translateTo) :
        pivotIndex(pivotIndex), subUnitIndex(subUnitIndex), translateTo(translateTo), backedUpPivotIndex(pivotIndex),
        backedUpTranslateToIndex(translateTo) {}

/**
 * @param subUnitIndex
 * @param pivotIndex   bead (universe index) placed at translateTo
 * @param translateTo  bead (universe index) whose position the pivot is moved to
 * @param indices      at least workingLimit bead indices
 * @param rng          generator for the random starting orientation
 */
SubUnit::SubUnit(unsigned int subUnitIndex, unsigned int pivotIndex, unsigned int translateTo, unsigned int workingLimit,
                 std::vector<unsigned int> &indices, PointSetModel *pModel, std::mt19937 &rng) :
        SubUnit(pivotIndex, subUnitIndex, translateTo) {

    if (indices.size() < workingLimit)
        throw std::invalid_argument("SubUnit: " + std::to_string(indices.size()) + " indices for a working limit of " + std::to_string(workingLimit));

    this->workingLimit = workingLimit;
    ensureCapacity(std::max(2u*workingLimit, 16u));
    std::copy(indices.begin(), indices.begin() + workingLimit, indicesInUse.begin());

    pPivotPoint = pModel->getBead(this->pivotIndex);
    pBackUpPivotPoint = pPivotPoint;
    pivotVector = pPivotPoint->getVec();
    translateToVector = pModel->getBead(this->translateTo)->getVec();
    backedUpCenteringVector = pivotVector;
    backedUpTranslationVector = translateToVector;

    for(unsigned int i=0; i < this->workingLimit; i++) coordinates[i] = pModel->getBead(indicesInUse[i])->getVec();

    // uniform random rotation: random unit quaternion (Shoemake), converted to a matrix
    std::uniform_real_distribution<float> u01(0.0f, 1.0f);
    const float u1 = u01(rng), u2 = u01(rng), u3 = u01(rng);
    const float a = std::sqrt(1.0f - u1), b = std::sqrt(u1);
    const float qw = a*std::sin(2.0f*(float)M_PI*u2), qx = a*std::cos(2.0f*(float)M_PI*u2);
    const float qy = b*std::sin(2.0f*(float)M_PI*u3), qz = b*std::cos(2.0f*(float)M_PI*u3);
    const float m[9] = {
            1-2*(qy*qy+qz*qz), 2*(qx*qy-qz*qw),   2*(qx*qz+qy*qw),
            2*(qx*qy+qz*qw),   1-2*(qx*qx+qz*qz), 2*(qy*qz-qx*qw),
            2*(qx*qz-qy*qw),   2*(qy*qz+qx*qw),   1-2*(qx*qx+qy*qy)};
    setRotationMatrix(m);
    transformAll();
    makeBackUp();
}

// row-major m (new = m * base) -> x1..z3 (new.x = x1 b.x + y1 b.y + z1 b.z, new.y = x2 ..., new.z = x3 ...)
void SubUnit::setRotationMatrix(const float m[9]) {
    x1 = m[0]; y1 = m[1]; z1 = m[2];
    x2 = m[3]; y2 = m[4]; z2 = m[5];
    x3 = m[6]; y3 = m[7]; z3 = m[8];
}

void SubUnit::getRotationMatrix(float m[9]) const {
    m[0] = x1; m[1] = y1; m[2] = z1;
    m[3] = x2; m[4] = y2; m[5] = z2;
    m[6] = x3; m[7] = y3; m[8] = z3;
}

void SubUnit::ensureCapacity(unsigned int needed) {
    if (needed <= capacity) return;
    unsigned int newCapacity = std::max(needed, std::max(2u*capacity, capacity + 16u));
    capacity = newCapacity;
    indicesInUse.resize(capacity);
    coordinates.resize(capacity);
    transformedCoordinates.resize(capacity);
    backedUpTransformedCoordinates.resize(capacity);
}

unsigned int SubUnit::locate(unsigned int beadIndex, const char * caller) const {
    auto end = indicesInUse.begin() + workingLimit;
    auto it = std::find(indicesInUse.begin(), end, beadIndex);
    if (it == end)
        throw std::logic_error(std::string("SubUnit::") + caller + ": bead " + std::to_string(beadIndex) +
                               " is not in subunit " + std::to_string(subUnitIndex) + " (" + std::to_string(workingLimit) + " beads in use)");
    return (unsigned int)std::distance(indicesInUse.begin(), it);
}

void SubUnit::makeBackUp() {
    backedUpPivotIndex = pivotIndex;
    backedUpTranslateToIndex = translateTo;
    pBackUpPivotPoint = pPivotPoint;

    backedUpCenteringVector = pivotVector;
    backedUpTranslationVector = translateToVector;

    x1_old = x1; x2_old = x2; x3_old = x3;
    y1_old = y1; y2_old = y2; y3_old = y3;
    z1_old = z1; z2_old = z2; z3_old = z3;

    std::copy(transformedCoordinates.begin(), transformedCoordinates.begin() + workingLimit, backedUpTransformedCoordinates.begin());
}

void SubUnit::transformPoint(const vector3 *base, vector3 *newVec) const {
    const vector3 b = *base - pivotVector;
    newVec->x = x1*b.x + y1*b.y + z1*b.z + translateToVector.x;
    newVec->y = x2*b.x + y2*b.y + z2*b.z + translateToVector.y;
    newVec->z = x3*b.x + y3*b.y + z3*b.z + translateToVector.z;
}

vector3 SubUnit::translateAndRotatePoint(const vector3 * originalVec) const {
    vector3 out;
    transformPoint(originalVec, &out);
    return out;
}

void SubUnit::transformAll() {
    for(unsigned int i=0; i<workingLimit; i++) transformPoint(&coordinates[i], &transformedCoordinates[i]);
}

/**
 * Undoing an add means decrementing workingLimit (the new bead is always last)
 */
void SubUnit::addIndex(unsigned int indexToAdd, PointSetModel *pModel){
    ensureCapacity(workingLimit + 1);
    indicesInUse[workingLimit] = indexToAdd;
    coordinates[workingLimit] = pModel->getBead(indexToAdd)->getVec();
    transformPoint(&coordinates[workingLimit], &transformedCoordinates[workingLimit]);
    ++workingLimit;
}

void SubUnit::undoAdd(){
    if (workingLimit == 0) throw std::logic_error("SubUnit::undoAdd: no bead to remove");
    workingLimit--;
}

/**
 * The removed bead is swapped to position workingLimit (just past the working set), so undoRemove restores it
 */
void SubUnit::removeIndex(unsigned int indexToRemove){
    const unsigned int pos = locate(indexToRemove, "removeIndex");
    workingLimit--;
    std::swap(indicesInUse[pos], indicesInUse[workingLimit]);
    std::swap(coordinates[pos], coordinates[workingLimit]);
    std::swap(transformedCoordinates[pos], transformedCoordinates[workingLimit]);
}

void SubUnit::undoRemove(){
    if (workingLimit >= capacity) throw std::logic_error("SubUnit::undoRemove: nothing was removed");
    workingLimit++;
}

/**
 * pivot -> universe position of bead translatedTo (untransformed)
 */
void SubUnit::translateToTransformedPosition(unsigned int translatedTo, PointSetModel * pModel){
    makeBackUp();
    translateTo = translatedTo;
    translateToVector = pModel->getBead(translatedTo)->getVec();
    transformAll();
}

void SubUnit::translateToPosition(const vector3 & position){
    makeBackUp();
    translateTo = UINT_MAX;
    translateToVector = position;
    transformAll();
}

void SubUnit::undoTranslateTo(PointSetModel *pModel) {
    pivotIndex = backedUpPivotIndex;
    translateTo = backedUpTranslateToIndex;
    pivotVector = backedUpCenteringVector;
    translateToVector = backedUpTranslationVector;
    pPivotPoint = (pModel != nullptr) ? pModel->getBead(backedUpPivotIndex) : pBackUpPivotPoint;
    std::copy(backedUpTransformedCoordinates.begin(), backedUpTransformedCoordinates.begin() + workingLimit, transformedCoordinates.begin());
}

/*
 * new pivot must be a bead of this subunit (checked); angles are small steps in degrees
 */
void SubUnit::rotateAt(unsigned int newPivotPoint, float angleInDegreesX, float angleInDegreesY, float angleInDegreesZ,
                       PointSetModel * pModel, bool inPlace){
    const unsigned int pos = locate(newPivotPoint, "rotateAt");
    makeBackUp();

    if (inPlace) translateToVector = transformedCoordinates[pos];   // the new pivot stays where it is now

    pivotIndex = newPivotPoint;
    pPivotPoint = pModel->getBead(newPivotPoint);
    pivotVector = pPivotPoint->getVec();

    // dR = Rz(c) Ry(b) Rx(a), then R <- dR R
    const float a = angleInDegreesX*convert, b = angleInDegreesY*convert, c = angleInDegreesZ*convert;
    const float ca = std::cos(a), sa = std::sin(a), cb = std::cos(b), sb = std::sin(b), cc = std::cos(c), sc = std::sin(c);
    const float d[9] = {
            cc*cb, cc*sb*sa - sc*ca, cc*sb*ca + sc*sa,
            sc*cb, sc*sb*sa + cc*ca, sc*sb*ca - cc*sa,
            -sb,   cb*sa,            cb*ca};
    float r[9], n[9];
    getRotationMatrix(r);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) n[3*i+j] = d[3*i]*r[j] + d[3*i+1]*r[3+j] + d[3*i+2]*r[6+j];
    // re-orthonormalise the rows (Gram-Schmidt) so that float round-off cannot accumulate over long runs
    float * r0 = n, * r1 = n + 3, * r2 = n + 6;
    float l = std::sqrt(r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2]); for(int k=0;k<3;k++) r0[k] /= l;
    float dp = r0[0]*r1[0] + r0[1]*r1[1] + r0[2]*r1[2]; for(int k=0;k<3;k++) r1[k] -= dp*r0[k];
    l = std::sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]); for(int k=0;k<3;k++) r1[k] /= l;
    r2[0] = r0[1]*r1[2] - r0[2]*r1[1]; r2[1] = r0[2]*r1[0] - r0[0]*r1[2]; r2[2] = r0[0]*r1[1] - r0[1]*r1[0];
    setRotationMatrix(n);
    transformAll();
}

void SubUnit::undoRotateTo() {
    x1 = x1_old; x2 = x2_old; x3 = x3_old;
    y1 = y1_old; y2 = y2_old; y3 = y3_old;
    z1 = z1_old; z2 = z2_old; z3 = z3_old;

    pivotVector = backedUpCenteringVector;
    translateToVector = backedUpTranslationVector;
    pPivotPoint = pBackUpPivotPoint;
    pivotIndex = backedUpPivotIndex;

    std::copy(backedUpTransformedCoordinates.begin(), backedUpTransformedCoordinates.begin() + workingLimit, transformedCoordinates.begin());
}

/*
 * replace bead indexOfCurrent by indexOfNew; the new bead is moved to the end of the working set (AddToPrSymX assumes
 * the newly added position is last)
 */
void SubUnit::swapIndices(unsigned int indexOfCurrent, unsigned int indexOfNew, PointSetModel * pModel) {
    swapLocale = locate(indexOfCurrent, "swapIndices");

    indicesInUse[swapLocale] = indexOfNew;
    coordinates[swapLocale] = pModel->getBead(indexOfNew)->getVec();
    transformPoint(&coordinates[swapLocale], &transformedCoordinates[swapLocale]);

    const unsigned int finalLocale = workingLimit - 1;
    std::swap(transformedCoordinates[swapLocale], transformedCoordinates[finalLocale]);
    std::swap(coordinates[swapLocale], coordinates[finalLocale]);
    std::swap(indicesInUse[swapLocale], indicesInUse[finalLocale]);
}
