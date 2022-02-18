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

#include "Anneal.h"

Anneal::Anneal(float highT,
               float percent,
               unsigned int highTempRounds,
               std::string fileprefix,
               float alpha,
               float beta,
               float eta,
               float lambda,
               float mu,
               unsigned int multiple,
               float accRate,
               float intercon) : highT(highT),
                                 percentAddRemove(percent),
                                 highTempRounds(highTempRounds),
                                 filenameprefix(std::move(fileprefix)),
                                 alpha(alpha),
                                 beta(beta),
                                 eta(eta),
                                 lambda(lambda),
                                 mu(mu),
                                 ccmultiple(multiple),
                                 asaAcceptanceRate (accRate),
                                 interconnectivityCutOff(intercon)
{


    //this->highTempStartForCooling = 0.00001; //0.00001
    this->highTempStartForCooling = 0.000016; //0.00001
    complementASAAcceptanceRate = 1.0f - accRate;
    intASAAcceptanceRate = (int)1000*accRate;
    intComplementASAAcceptanceRate = (int)1000*complementASAAcceptanceRate;

    contactsDistribution.resize(13);
    /*
     * this distribution is set to q-max 0.4 and bead_radius = 0.5*bin_width
     */
    distributionlimit=13;
    contactsDistribution[0] = 0.0;
//    contactsDistribution[1] = 0.25419d;
    contactsDistribution[1] = 0.27;
    contactsDistribution[2] = 0.43;
//    contactsDistribution[2] = 0.317564d;
    contactsDistribution[3] = 0.241696;
    contactsDistribution[4] = 0.117231;
    contactsDistribution[5] = 0.0447598;
    contactsDistribution[6] = 0.0190639;
    contactsDistribution[7] = 0.00549451;
    contactsDistribution[8] = 0.0;
    contactsDistribution[9] = 0.0;
    contactsDistribution[10] = 0.0;
    contactsDistribution[11] = 0.0;
    contactsDistribution[12] = 0.0;

    double totaltemp = 0;
    for(unsigned int i=0; i<distributionlimit; i++){
        totaltemp += contactsDistribution[i];
    }

    for(unsigned int i=0; i<distributionlimit; i++){ // normalize
        contactsDistribution[i] *= 1.0/totaltemp;
    }
}

/**
 *
 * @param flags
 * @param bead_indices
 * @param upTo is the workingLimit of beads_indices
 * @param pModel
 * @return
 */
float Anneal::calculateCVXHULLVolume(char *flags, std::vector<unsigned int> *bead_indices, const unsigned int upTo, PointSetModel *pModel) {

    unsigned int numpoints = 3*upTo;
    coordT points[numpoints];

    const unsigned int * ptr = (*bead_indices).data();

    for (unsigned int i=0; i<upTo; i++){
        beadToPoint(&points[i*3], pModel->getBead(ptr[i]));
    }

    qh_new_qhull (3, upTo, points, 0, flags, nullptr, nullptr);
    float volume_test = (float)(qh totvol);

    //qh totarea;
    qh_freeqhull(true);

    return volume_test;
}

float Anneal::surfaceToVolume(const unsigned int total, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates){

    POWERSASA::PowerSasa<float, Eigen::Vector3f> *ps =
            new POWERSASA::PowerSasa<float, Eigen::Vector3f>(coordinates, weights, 1, 1, 1, 1);

    ps->calc_sasa_all();
    float volume = 0.0;
    float sasa = 0.0;
    for (unsigned int i = 0; i < total; ++i) {
        sasa += ps->getSasa()[i];
        volume += ps->getVol()[i];
    }

    delete ps; //maximizes surface area while minimizing number of regions
    return (float)(36*M_PI*(volume*volume)/(sasa*sasa*sasa)); //sphericity : minimizing moves away from a sphere
    // lots of surface area supports an extended model
    // minimal surface area supports a compact object
    //return sasa;
};


/**
 * Creates a distribution of contacts from selected lattice points specified by beads_in_use
 *
 * @param distribution
 * @param beads_in_use
 * @param pModel
 */
void Anneal::populateContactsDistribution(std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, PointSetModel *pModel){

    for(int i=0; i<13; i++){
        distribution[i]=0.0d;
    }

    for (auto const & it : *beads_in_use){
        ++distribution[ numberOfContactsFromSet(beads_in_use, pModel, it) ];
    }
}

unsigned int Anneal::numberOfContactsFromSet(std::set<unsigned int> *beads_in_use,
                                             PointSetModel *pModel,
                                             unsigned int const selectedIndex){

    auto it = pModel->getPointerToNeighborhood(selectedIndex);
    unsigned int neighborContacts = 0;

    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighhborLimit = pModel->getNeighborLimit();

    for (unsigned int i=0; i< totalNeighbors; i++){

        unsigned int neighbor = *(it+i);

        if ((neighbor < neighhborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){
            neighborContacts += 1;
        } else if (neighbor == neighhborLimit) {
            break;
        }
    }

    return neighborContacts;
}

/**
 * pre-calculation, populate vector of bins for all pairwise distances
 *
 * @param pModel
 * @param pData
 */
void Anneal::fillPrBinsAndAssignTotalBin(PointSetModel * pModel, PofRData * pData){

    maxbin = pModel->populateBins(pData); // populate the entire over entire Universe

    maxbin += 1;

    totalBins = pData->getShannonBins(); // set by experimental input Data file
    // maxbin and totalBins will be different
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

}


/**
 * Calcualte P(r) distribution using selected set of bead_indices (bead indices must be sorted!)
 * binCount is a small vector
 * Should not use if max index of pModelBin is greater than unsigned_int max
 * @param bead_indices
 * @param binCount
 * @param workingLimit
 * @param totalBeadsInSphere
 * @param pModel
 * @param pData
 * @return
 */
void Anneal::calculateModelPrDistributionDirect(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, const unsigned int workingLimit, PointSetModel *pModel, PofRData *pData){

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);
    // calculate P(r) for beads
    // calculate distances and determine bin
    const unsigned int * const pBead = (*bead_indices).data();

    for(unsigned int m=0; m < workingLimit; m++){ // set row

        const vector3 * firstVec = &(pModel->getBead(pBead[m])->getVec());

        for(unsigned int n=(m+1); n < workingLimit; n++){ // iterate over columns
            (*binCount)[ pData->convertToBin((*firstVec - pModel->getBead(pBead[n])->getVec()).length()) ]++;
        }
    }
}

/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSymCE(std::vector<unsigned int> *subUnit_indices,
                                               std::set<unsigned int> & subUnit_indices_tree,
                                               std::map<unsigned int, std::vector<vector3> > &map,
                                               std::vector<double> & contactsDistributionOfModel,
                                               std::vector<unsigned int> *binCount,
                                               const unsigned int indicesWorkingLimit,
                                               unsigned int &violations,
                                               PointSetModel *pModel, PofRData *pData) {

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3 *> coordinates(totalCoordinates);

    vector3 ** const pCoords = coordinates.data();
    unsigned int * const pSubInd = (*subUnit_indices).data();
    // create first subunit from selected indices and populate coordinates
//    for (unsigned int i=0; i < indicesWorkingLimit; i++){
//        tempBead = pModel->getBead(pSubInd[i]);
//        ptr[i] = vector3(tempBead->getX(), tempBead->getY(), tempBead->getZ());
//    }

    unsigned int index=0;
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        auto mit = map.find(pSubInd[i]);
        vector3 * const pVecBase = (*mit).second.data();

        for(unsigned int s=0; s<pModel->getNumberOfSubUnits(); s++){ // iterate through map, add all symmetry related points
            pCoords[index] = &pVecBase[s];
            index++;
        }
    }

//    unsigned int count = indicesWorkingLimit;
//    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
//        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
//    }

    // calculate Pr, order is irrelavant (all pairwise)
    // parallelize this code

    float distance_to;
    for (unsigned int i=0; i < totalCoordinates; i++) {

        const vector3 * tempVec1 = pCoords[i];
        /*
         * each thread to have its own copy of tempVec1
         */

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = (*tempVec1 - *pCoords[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // unrolling could lead to a RACE condition of updating same location at same time
        }
    }


    /*
     * for each lattice point within parent, determine how many violations it makes with other points in symmetry mates
     */
//    violation=0;
//    for (unsigned int i=0; i < indicesWorkingLimit; i++) {
//
//        const vector3 * ptempVec = &coordinates[i];
//
//        for (unsigned int next_i = indicesWorkingLimit; next_i < totalCoordinates; next_i++) {
//            // calculate distance and convert to bin
//            distance_to = ((*ptempVec) - coordinates[next_i]).length();
//            if (distance_to < violation_limit){
//                violation++;
//            }
//        }
//    }

    std::fill(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), 0.0d);
    violations=0;
    unsigned int tempContacts;
    for(unsigned int ind=0; ind<indicesWorkingLimit; ind++){
        tempContacts = numberOfContactsFromSet(&subUnit_indices_tree, pModel, pSubInd[ind]);
        const vector3 *baseVec = pCoords[ind];

        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = (*baseVec - *pCoords[s]).length();
            if (distance_to < contactCutOff){ // since violation_limit is always < contactCutOff, only check if true
                tempContacts++;
                if (distance_to < violation_limit){
                    violations++;
                }
            }
        }
        tempContacts = (tempContacts > 12) ? 12 : tempContacts;
        ++contactsDistributionOfModel[tempContacts];
    }
}

/**
 * bead_indices should be sorted to make insertion into hull as efficient as possible
 * finds all unused points within CVX hull including neighbors of hull vertices
 * @param workingLimit
 * @param bead_indices
 * @param hull
 * @param pModel
 * @param totalBeadsInSphere
 */
void Anneal::recalculateDeadLimit(unsigned int workingLimit, std::vector<unsigned int> &bead_indices, std::set<unsigned int> &hull, PointSetModel * pModel){

    hull.clear();

    pointT testPoint[3];
    boolT isoutside;
    realT bestdist;
    char flags[] = "qhull FA";

    coordT hullPoints2[3*workingLimit];
    // calculate CVX Hull from selected indices
    const auto num = (unsigned int)bead_indices.size();
    const unsigned int * ptr = (num > 0) ? bead_indices.data() : nullptr;
    std::vector<unsigned int> active_indices(workingLimit);

    for (unsigned int i = 0; i < workingLimit; i++) {
        beadToPoint(&hullPoints2[i*3], pModel->getBead(ptr[i]));
        active_indices[i] =ptr[i];
    }

    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, nullptr, nullptr);

    for(unsigned int i = workingLimit; i < num; i++) {
        unsigned int value = ptr[i];
        beadToPoint(testPoint, pModel->getBead(value));
        // exclude HULL points, for each bead, determine if outside HULL
        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);

        if (!isoutside){
            hull.insert(value); // if not outside, keep it
        }
    }

    std::set<unsigned int> beads_in_use(bead_indices.begin(), bead_indices.begin() + workingLimit);

    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    unsigned int neighbor;
    vertexT * vertices = qh vertex_list;
    auto totalV = (unsigned int)qh num_vertices;


    for (unsigned int v = 0; v < totalV; v++) { // check vertices of CVX Hull

        auto pptr = pModel->getDirectPointerToNeighborhood();
        unsigned int location = totalNeighbors*active_indices[qh_pointid( vertices->point)];

        for (unsigned int j=0; j < totalNeighbors; j++){ // add neighbors if not in use
            neighbor = pptr[location+j];
            if (neighbor < neighborLimit && beads_in_use.find(neighbor) == beads_in_use.end()){
                hull.insert(neighbor);
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
        vertices = vertices->next;
    }

    qh_freeqhull(true);
}

/**
 * Calcualte P(r) distribution using selected set of bead_indices (bead indices must be sorted!)
 * binCount is a small vector
 * Should not use if max index of pModelBin is greater than unsigned_int max
 * @param bead_indices
 * @param binCount
 * @param workingLimit
 * @param totalBeadsInSphere
 * @param pModel
 * @return
 */
void Anneal::calculateModelPrDistribution(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, unsigned int workingLimit, unsigned int totalBeadsInSphere, PointSetModel *pModel) {

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);

    unsigned int row;
    unsigned int row2;
    const unsigned short int * pModelBin = pModel->getPointerToBins();

    //const int * pBeadIndices = &(*bead_indices)[0];
    const unsigned int * ptr = (*bead_indices).data();

    // calculate P(r) for beads
    for(unsigned int m=0; m < workingLimit; m++){ // set row

        //row = (unsigned int)*(pBeadIndices + m);
        row = ptr[m];
        /*
         * if row is zero, row2 can not be less than zero.
         * move -1 into the loop, insures number is never assigned non-negative
         * Underflow error possible
         */
        row2 = (row*totalBeadsInSphere) - (row*(row+1))/2 - row;// - 1;

        for(unsigned int n=(m+1); n < workingLimit; n++){ // iterate over columns
            (*binCount)[ *( pModelBin + ((row2 + ptr[n]) - 1) ) ]++;
        }
    }

}


/**
 * for each selected lattice position within workingLimit
 * grab lattice points that comprise its neighborhood
 * and for each point not already within workingLimit, move to within hull sett
 */
void Anneal::populateLayeredDeadlimitUsingSet(std::set<unsigned int> & beads_in_use, std::set<unsigned int> & hull, PointSetModel * pModel) {

    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    unsigned int neighbor;
    hull.clear();

    for(const auto & it : beads_in_use){
        //hull.insert(it);
//        auto neighborhood = pModel->getPointerToNeighborhood(it);
        const auto ptr = pModel->getDirectPointerToNeighborhood();
        unsigned int location = totalNeighbors*it;

        for (unsigned int j=0; j < totalNeighbors; j++){
            neighbor = ptr[location+j];

            if (neighbor < neighborLimit && beads_in_use.find(neighbor) == beads_in_use.end()){
                hull.insert(neighbor);
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
    }
}

/**
 * from Vincent A. Cicirello
 * On the Design of an Adaptive Simulated Annealing Algorithm
 *
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void Anneal::updateASAConstantTemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp){


    double stepEval = (double)index/evalMax;
    double lamRate=asaAcceptanceRate;

    if (stepEval < 0.15) {
        //lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
        lamRate = asaAcceptanceRate+complementASAAcceptanceRate*pow(intComplementASAAcceptanceRate, -stepEval*6.666667);
    } else if (stepEval >= 0.15){
        lamRate = asaAcceptanceRate;
    }

    if (acceptRate > lamRate){
        temp = 0.999*temp;
    } else {
        temp = temp*1.001001001001;
    }

    inv_temp = 1.0/temp;
}

/**
 * from Vincent A. Cicirello
 * On the Design of an Adaptive Simulated Annealing Algorithm
 *
 * @param index
 * @param evalMax
 * @param acceptRate
 * @param temp
 * @param inv_temp
 */
void Anneal::updateASATemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp){

    double stepEval = (double)index/evalMax;
    double lamRate=asaAcceptanceRate;

    if (stepEval < 0.15) {
        //lamRate = 0.44+0.56*pow(560, -stepEval*6.666667);
        lamRate = asaAcceptanceRate+complementASAAcceptanceRate*pow(intComplementASAAcceptanceRate, -stepEval*6.666667);
    } else if (stepEval >= 0.15 && stepEval < 0.65){
        lamRate = asaAcceptanceRate;
    } else if (stepEval >= 0.65){
        lamRate = asaAcceptanceRate*pow(intASAAcceptanceRate, -(stepEval - 0.65)*2.857142857);
    }

    temp = (acceptRate > lamRate) ? (0.999*temp) : (temp*1.001001001001);

    inv_temp = 1.0/temp;
}

// returns points of the convex hull within the hullpts set
void Anneal::getHullPoints(std::set<unsigned int> &hullpts, std::set<unsigned int> &beads_in_use, PointSetModel * pModel){
    unsigned int workingLimit = beads_in_use.size();
    char flags[] = "qhull FA";
    hullpts.clear();

    coordT hullPoints2[3*workingLimit];
    std::vector<unsigned int> active_indices(workingLimit);
    // calculate CVX Hull from selected indices
    unsigned int count=0;
    for(const auto & ind : beads_in_use){
        beadToPoint(&hullPoints2[count*3], pModel->getBead(ind));
        active_indices[count] =ind;
        count++;
    }
    // calculate convex hull
    qh_new_qhull(3, workingLimit, hullPoints2, 0, flags, nullptr, nullptr);

    vertexT * vertices = qh vertex_list;
    auto totalV = (unsigned int)qh num_vertices;

    // only move CVX hull points
    for (unsigned int v = 0; v < totalV; v++) { //
        hullpts.insert(active_indices[qh_pointid( vertices->point)]);
        vertices = vertices->next;
    }

    qh_freeqhull(true);
}


/*!
 *
 */
double Anneal::connectivityPotential(unsigned int numberOfComponents){

    switch(numberOfComponents) {
        case 1:
            return 0.0d;
        case 2:
            return 0.1d;
        case 3:
            return 0.30d;
        case 4:
            return 0.70d;
        case 5:
            return 100000.0d;
        case 6:
            return 1000000.0d;
        default:
            return 1000000.0d*numberOfComponents;
    }
}

/**
 * beads_in_use must already contain addMe
 * @param addMe
 * @param distribution
 * @param beads_in_use
 * @param pModel
 */
void Anneal::addToContactsDistribution(unsigned int addMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, PointSetModel *pModel){
    // addMe will have n-number of contacts and for each
    // get all neighbors for addMe
    double * const pDistribution = distribution.data(); // initialized as empty in Model class
    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(addMe);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // now adjust counts
    for(unsigned int i=0; i<count; i++){
        unsigned int tempNumber=numberOfContactsFromSet(beads_in_use, pModel, neighborsInuse[i]);
        --pDistribution[ tempNumber-1 ];
        ++pDistribution[ tempNumber ];
//        --distribution[ tempNumber-1 ]; // remove prior contribution
//        ++distribution[ tempNumber ];    // add new contribution
    }
    // now do addMe
    ++pDistribution[numberOfContactsFromSet(beads_in_use, pModel, addMe)];
}


/**
 * beads_in_use must already contain removeMe
 * @param addMe
 * @param distribution
 * @param beads_in_use
 * @param pModel
 */
void Anneal::removeFromContactsDistribution(unsigned int removeMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, PointSetModel *pModel){
    // get all neighbors for removeMe
    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(removeMe);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // now adjust counts

    for(unsigned int i=0; i<count; i++){
        unsigned int currentNumber=numberOfContactsFromSet(beads_in_use, pModel, neighborsInuse[i]);
        --distribution[ currentNumber ];      // remove prior contribution
        ++distribution[ currentNumber-1 ];    // add new contribution
    }

    // now remove contributions from removeMe bead
    --distribution[ numberOfContactsFromSet(beads_in_use, pModel, removeMe)];
}



/**
 * Given the selectedIndex, compile list of indices that are direct neighbors in use
 * Pick one at random
 *
 */
unsigned int Anneal::getConnectedNeighborFromSet(std::set<unsigned int> *beads_in_use,
                                                 PointSetModel *pModel,
                                                 unsigned int & selectedIndex){

    const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> neighborsInuse(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = *(it+i);
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) != beads_in_use->end()){ // if end of set, means not in use
            neighborsInuse[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return neighborLimit;
    } else if (count == 1){
        return neighborsInuse[0];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);

    return neighborsInuse[randomIndex(gen)];
}


/**
 * for selected Index, calculate a probability based on number of contacts
 * if bead has few contacts, it will have lots of available neighbors
 *
 */
float Anneal::getProbabilityOfSelectedIndex(std::set<unsigned int> *beads_in_use,
                                                PointSetModel *pModel,
                                                unsigned int & selectedIndex){

    const unsigned int * ptr = pModel->getDirectPointerToNeighborhood();
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();

    unsigned int countAvailableNeighbors=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        //unsigned int neighbor = *(it+i);
        unsigned int neighbor = ptr[totalNeighbors*selectedIndex + i];
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) == beads_in_use->end()){ // if end of set, means not in use
            countAvailableNeighbors++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    return (float)countAvailableNeighbors/12.0f;
}



/**
 * pick random point that is in use
 * get set of neighbors that are available
 * randomly select one
 */
unsigned int Anneal::getUseableNeighborWeighted(std::set<unsigned int> *beads_in_use,
                                               PointSetModel *pModel,
                                               unsigned int & selectedIndex){


    const unsigned int * ptr = pModel->getDirectPointerToNeighborhood();
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    unsigned int countAvailable=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = ptr[totalNeighbors*selectedIndex + i];
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) == beads_in_use->end()){ // if end of set, means not in use
            possibleNeighbors[countAvailable] = neighbor;
            countAvailable++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> randomProb(0,1);

    if (randomProb(gen) < (double)countAvailable/((double)totalNeighbors)){
        std::uniform_int_distribution<unsigned int> randomIndex(0,countAvailable-1);
        return possibleNeighbors[randomIndex(gen)];
    } else {
        return neighborLimit;
    }
}

/**
 * pick random point that is in use
 * get set of neighbors that are available
 * randomly select one
 */
unsigned int Anneal::getUseableNeighborFromSetCE(std::set<unsigned int> *beads_in_use,
                                               PointSetModel *pModel,
                                               unsigned int & selectedIndex){


    const unsigned int * ptr = pModel->getDirectPointerToNeighborhood();
    //const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        unsigned int neighbor = ptr[totalNeighbors*selectedIndex + i];
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) == beads_in_use->end()){ // if end of set, means not in use
            possibleNeighbors[count] = neighbor;
            count++;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return neighborLimit;
    } else if (count == 1){
        return possibleNeighbors[0];
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);

    return possibleNeighbors[randomIndex(gen)];
}


/**
 * pick random point that is in use
 * get set of neighbors that are available
 * randomly select one
 */
unsigned int Anneal::getUseableNeighborFromSet(std::set<unsigned int> *beads_in_use,
                                               PointSetModel *pModel,
                                               unsigned int & selectedIndex){


    const unsigned int * ptr = pModel->getDirectPointerToNeighborhood();
    //const std::vector<unsigned int>::iterator it = pModel->getPointerToNeighborhood(selectedIndex);
    // go through each member of the neighborhood
    // determine their current energy state and after if bead is moved
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    std::vector<unsigned int> possibleNeighbors(totalNeighbors);

    unsigned int count=0;

    for (unsigned int i=0; i< totalNeighbors; i++){
        //unsigned int neighbor = *(it+i);
        unsigned int neighbor = ptr[totalNeighbors*selectedIndex + i];
        if ((neighbor < neighborLimit ) && beads_in_use->find(neighbor) == beads_in_use->end()){ // if end of set, means not in use
            possibleNeighbors[count] = neighbor;
            count++;
        } else if (neighbor == neighborLimit) {
            break;
        }
    }

    // what happens if no neighbors?
    if (count==0){
        return neighborLimit;
    } else if (count == 1){
        return possibleNeighbors[0];
    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<unsigned int> randomIndex(0,count-1);

    return possibleNeighbors[randomIndex(gen)];
}


/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSym(std::vector<unsigned int> *subUnit_indices, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, unsigned int &violation,  PointSetModel *pModel, PofRData *pData) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    vector3 * const ptr = coordinates.data();
    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        ptr[i] = pModel->getBead((*subUnit_indices)[i])->getVec();
    }


    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate Pr, order is irrelavant
    float distance_to;
    violation=0;
    for (unsigned int i=0; i < indicesWorkingLimit; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < indicesWorkingLimit; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = indicesWorkingLimit; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
            if (distance_to < violation_limit){
                violation++;
            }
        }
    }

    for (unsigned int i=indicesWorkingLimit; i < totalCoordinates; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }
}


/**
 * beadsInUse does not have to be sorted
 *
 * prBins is the model P(r) distribution
 *
 * returns the number of violations adjusted with respect to removeMeSubUnitIndex
 *
 *
 */
unsigned int Anneal::removeFromPrSym(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, PointSetModel *pModel, PofRData *pData){

    unsigned int findIt=0;
    // find index of bead to remove
    const unsigned int * const ptr = beadsInUse.data();
    for(unsigned int i=0; i < workingLimit; i++){
        if (ptr[i] == removeMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const unsigned int stopAt = findIt;
    //
    // remove interdomain distances
    //
    float distance_to;

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < workingLimit; i++){
        coordinates[i] = pModel->getBead(ptr[i])->getVec();
    }

    unsigned int count = workingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }


    // adjust Pr
    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }


    /**
     *
     * we double remove(count) the symmetry related vectors, need to add back
     */
    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int ss=s+1; ss < subUnits; ss++){
            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

    return 0;
}

/**
* The new position to Add MUST BE WITHIN THE WORKING LIMIT
* beadsInUse must be sorted up to workingLimit
*
* @param addMeSubUnitIndex
* @param beadsInUse
* @param workingLimit
* @param prBins
* @param pModel
* @param pData
* @return
*/
unsigned int Anneal::addToPrSym(unsigned int addMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, PointSetModel *pModel, PofRData *pData){
    // each bead has a sym related partner in pModel
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin

    unsigned int findIt=0;
    const unsigned int * const ptr = beadsInUse.data();
    for(unsigned int i=(workingLimit-1); i>=0; i--){
        if (ptr[i] == addMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const int stopAt = findIt;
    // add interdomain distances
    //
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * const pC = coordinates.data();

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < workingLimit; i++){
        pC[i] = pModel->getBead(ptr[i])->getVec();
    }

    unsigned int count = workingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }

    // adjust Pr
    float distance_to;
    for (unsigned int s=0; s < subUnits; s++){  //
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &pC[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - pC[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - pC[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

    /*
     * correct for double counting
     */
    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &pC[basis];

        for (unsigned int ss=s+1; ss < subUnits; ss++){
            distance_to = ((*prefVec) - pC[stopAt+ss*workingLimit]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }
    /*
     * Should include contacts distribution
     *
     *
     * count any symmetry related clashes
     */
    unsigned int totalviolations=0;
    for (unsigned int i=0; i < workingLimit; i++){
        const vector3 * pBase = &pC[i];
        for (unsigned int next=workingLimit; next < totalCoordinates; next++){
            distance_to = ((*pBase) - pC[next]).length();
            if (distance_to < violation_limit){
                ++totalviolations;
            }
        }
    }

    return totalviolations;
}

/**
 * count the number of violations of parent subunit
 * @param bead_indices
 * @param beadIndiciesWorkingLimit
 * @param pModel
 * @param pData
 * @return
 */
unsigned int Anneal::getViolations(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, PointSetModel *pModel){

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

//    vector3 * tempVec1;
//    Bead * tempBead;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
//        tempBead = pModel->getBead((*subUnit_indices)[i]);
        coordinates[i] = pModel->getBead((*subUnit_indices)[i])->getVec(); // copy assignment
//        tempVec1 = &coordinates[i];
//        (*tempVec1).x = tempBead->getX();
//        (*tempVec1).y = tempBead->getY();
//        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

//    std::string nameOf = "violations.pdb";
//    FILE * pFile = fopen(nameOf.c_str(), "a");

    float distance_to;
    unsigned int violations = 0;
    // count violations
    for(unsigned int i=0; i<indicesWorkingLimit; i++) {

        const vector3 *baseVec = &coordinates[i];
        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = ((*baseVec) - coordinates[s]).length();
            if (distance_to < violation_limit){
//                pModel->printAtomLine(pFile, s, "A", s, coordinates[s].x, coordinates[s].y, coordinates[s].z);
                violations++;
            }
        }
    }

//    fclose(pFile);
    // check for connectivity
    return violations;
}

/**
 * Calculate both the number of violations (clashes) from symmetry mates and the contact distirbution ofthe subunit
 * @param subUnit_indices (sorted due to SET container)
 * @param contactsDistributionOfModel
 * @param indicesWorkingLimit
 * @param pModel
 * @return
 */
unsigned int Anneal::getViolationsFromSet(std::set<unsigned int> *subUnit_indices, std::vector<double> & contactsDistributionOfModel, const unsigned int indicesWorkingLimit, PointSetModel *pModel){


    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);


    // create first subunit from selected indices and populate coordinates
    vector3 * const ptr = coordinates.data();
    unsigned int index=0;
    for(auto & ind : *subUnit_indices){ // sorted in order because of set
        ptr[index] = pModel->getBead(ind)->getVec();
        index++;
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }


//    std::string nameOf = "violations.pdb";
//    FILE * pFile = fopen(nameOf.c_str(), "a");

    float distance_to;
    unsigned int violations = 0;
    std::fill(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), 0.0d);
    // count violations
    index=0;
    for(auto & ind : *subUnit_indices){
        unsigned int tempContacts = numberOfContactsFromSet(subUnit_indices, pModel, ind);
        const vector3 *baseVec = &ptr[index];

        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = ((*baseVec) - ptr[s]).length();

            //           if (distance_to < contactCutOff){ // this includes additional symmetry mate contacts
            //               tempContacts++;
            if (distance_to < violation_limit){
//                    pModel->printAtomLine(pFile, s, "A", s, ptr[s].x, ptr[s].y, ptr[s].z);
                violations++;
            }
            //           }
        }
        index++;
        //       tempContacts = (tempContacts > 12) ? 12 : tempContacts;
        ++contactsDistributionOfModel[tempContacts];
    }

//    fclose(pFile);
    // entire Symmetry made molecule
//    for(unsigned int index=0;index<totalCoordinates; index++){
//        const vector3 *baseVec = &ptr[index];
//        unsigned int counter=0;
//        for(unsigned int index2=index+1;index2<totalCoordinates; index2++){
//            const vector3 *baseVec2 = &ptr[index2];
//            distance_to = (*baseVec - *baseVec2).length();
//            if (distance_to < contactCutOff){
//                counter++;
//                if (counter >= 12){
//                    break;
//                }
//            }
//        }
//        ++contactsDistributionOfModel[counter];
//    }


//    for(auto & ind : *subUnit_indices){
//        unsigned int tempContacts = numberOfContactsFromSet(subUnit_indices, pModel, ind);
//        const vector3 *baseVec = &ptr[index];
//
//        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
//            distance_to = ((*baseVec) - ptr[s]).length();
//
//            if (distance_to < contactCutOff){
//                tempContacts++;
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//            }
//        }
//        index++;
//
//        tempContacts = (tempContacts > 12) ? 12 : tempContacts;
//        ++contactsDistributionOfModel[tempContacts];
//    }

    // check for connectivity
    return violations;
}



void Anneal::calculateModelParametersSymmetry(std::set<unsigned int> *subUnit_indices, PointSetModel *pModel) {

    //const unsigned int subUnits = pModel->getNumberOfSubUnits();

    const vector3 vecXY(1,1,0);
    const vector3 vecZ(0,0,1);

    float maxXY=0, maxZ=0;

    for(auto & ind : *subUnit_indices){ // sorted in order because of set
        const vector3 * tempVec1 = &(pModel->getBead(ind)->getVec());
        vector3 temp = (*tempVec1)*vecXY;
        if (temp.length() > maxXY){
            maxXY=temp.length();
        }

        temp = (*tempVec1)*vecZ;
        if (temp.length() > maxZ){
            maxZ=temp.length();
        }
    }

    logger("SEARCH SPACE PARAMETERS", "FOR A CYLINDER");
    logger("HEIGHT", formatNumber(2*(maxZ + pModel->getBeadRadius()),2));
    logger("RADIUS", formatNumber(maxXY + pModel->getBeadRadius(),2));
}


void Anneal::printParameters(std::vector<float> * accept, std::vector<double> * temp, std::vector<float> * divergence, std::vector<unsigned int> * wl){

    int total = accept->size();

    FILE * pFile;
    const char *outputFileName;
    std::string nameOf = "run_parameters.txt";
    outputFileName = nameOf.c_str() ;
    pFile = std::fopen(outputFileName, "w");

    std::string index;

    for (int i=0; i < total; i++){
        fprintf(pFile, "%i %.5E %0.9f %0.9f %i\n", (i+1), (*accept)[i], (*temp)[i], (*divergence)[i], (*wl)[i] );
    }

    fclose(pFile);
}

bool Anneal::setAnchorPoints(std::string anchorFileName, std::string pdbFile, PointSetModel *pModel){

    PDBModel pdbModel(pdbFile, true, false); // centered Coordinates

    // if anchor points are in pdbFile, return true, else return false
    // CHAIN, RESIDUE NUMBER, ATOM?
    // ATOM     54  O   GLY A   8
    const unsigned int totalAtoms = pdbModel.getTotalCoordinates();
    std::string line;

    std::ifstream anchorFile (anchorFileName.c_str());
    boost::regex pdbStart("ATOM");
    boost::regex residue("RESID");
    boost::regex lineFormat("\\w+\\s+[0-9]+\\s+\\w+[A-Z0-9]+", boost::regex::icase);
    boost::regex component_id("COMPONENT_ID");
    boost::regex volume("VOLUME");
    boost::regex chain("CHAIN");
    boost::regex wat("HOH");
    boost::regex hash("#");

    auto pdbResIDs = pdbModel.getResIDIterator();
    std::vector<std::string>::const_iterator pdbAtomTypes = pdbModel.getAtomTypeIterator();
    std::vector<std::string>::const_iterator pdbChainIds = pdbModel.getChainIDIterator();
    const unsigned int totalBeads = pModel->getTotalNumberOfBeadsInUniverse();

    Bead * currentBead;

    // find closest non-seed bead position!
    // format of Anchor file
    std::vector<std::string> tempLine;
    std::vector<std::string> splitLine;
    std::vector<unsigned int> resids;
    std::vector<float> volumes;
    std::vector<std::string> ids;

    std::string currentComponentID;

    // get lines in the file
    if (anchorFile.is_open()) {
        while(!anchorFile.eof()) {
            std::getline(anchorFile, line);
            boost::algorithm::trim(line);
            tempLine.push_back(line);
        }
    }
    anchorFile.close();

    /*
     * get COMPONENT_IDs and VOLUMEs
     * ex: COMPONENT_ID B VOLUME 8000
     */
    try {
        for(auto it = tempLine.begin(); it != tempLine.end(); ++it) {

            boost::split(splitLine, *it, boost::is_any_of("\t  "), boost::token_compress_on);

            if (!(*it).empty() && boost::regex_search(splitLine[0], component_id) && splitLine.size() == 4 && boost::regex_search(*it, volume)){
                components.emplace_back( Component(splitLine[1], stof(splitLine[3]), pModel) );
                volumes.push_back(stof(splitLine[3]));
            } else if ( boost::regex_search(splitLine[0], component_id) && !splitLine[0].empty()) {
                throw std::invalid_argument( "COMPONENT ID or VOLUME NOT SPECIFIED : \n\t" + *it  + " \n");
            }
        }

    } catch (std::exception &err) {
        std::cerr<<"Caught "<<err.what()<<std::endl;
        std::cerr<<"Type "<<typeid(err).name()<<std::endl;
        exit(0);
    }

    // for each component add the resids
    try {
        for(auto it = tempLine.begin(); it != tempLine.end(); ++it) {

            boost::split(splitLine, *it, boost::is_any_of("\t  "), boost::token_compress_on);

            if (!(*it).empty() && boost::regex_search(splitLine[0], residue) && splitLine.size() == 6 && boost::regex_search(*it, component_id) && boost::regex_search(*it, chain)){

                // component_ID must be in the list, if not throw exception
                std::string tempId = splitLine[5];
                auto fit = std::find_if(components.begin(), components.end(), [&tempId](const Component& obj) {return obj.getID() == tempId;});

                if (fit != components.end()) {
                    // found element. it is an iterator to the first matching element.
                    // if you really need the index, you can also get it:
                    //auto index = std::distance(components.begin(), fit);
                    auto tempResid = (unsigned int) std::stoi(splitLine[1]);
                    if (tempResid > 1){ // check if RESID is in PDB model
                        (*fit).addResid(tempResid, splitLine[3]);
                    } else {
                        throw std::invalid_argument( "IMPROPER RESID: \n\t" + *it  + " RESID \n" + std::to_string(tempResid) + " \n");
                    }
                } else {
                    throw std::invalid_argument( "COMPONENT ID MISSING OR INCORRECT: \n\t" + *it  + " \n");
                }

            } else if ( !(*it).empty() && boost::regex_search(splitLine[0], residue) && splitLine.size() < 6 && boost::regex_search(*it, component_id) ) {
                throw std::invalid_argument( "COMPONENT ID or RESID NOT SPECIFIED : \n\t" + *it  + " \n");
            }
        }
    } catch (std::exception &err) {
        std::cerr<<"Caught "<<err.what()<<std::endl;
        std::cerr<<"Type "<<typeid(err).name()<<std::endl;
        exit(0);
    }

    // for each Component, find lattice point that is central to the residue
    // map resid to structure, for each resid, grab all the atoms and calculate average
    float xpos, ypos, zpos;
    //float b2 = pModel->getBeadRadius()*pModel->getBeadRadius()*1.13;
    float diffx, diffy, diffz;

    for(auto it = components.begin(); it != components.end(); ++it) {

        float dis2;
        for(unsigned int r_i=0; r_i < (it->getTotalResids()); r_i++){
            xpos=0;
            ypos=0;
            zpos=0;
            unsigned int atomCounter=0;

            for (unsigned int i=0; i < totalAtoms; i++){ // calculate average position of residue
                // get all atom positions with matching resid and chain
                if ( it->getResidByIndex(r_i) == *(pdbResIDs + i) && (it->getChainsByIndex(r_i) == *(pdbChainIds + i)) ) {
                    xpos += *(pdbModel.getCenteredXVec() + i);
                    ypos += *(pdbModel.getCenteredYVec() + i);
                    zpos += *(pdbModel.getCenteredZVec() + i);
                    atomCounter++;
                } else if (it->getResidByIndex(r_i) < *(pdbResIDs + i)) {
                    break;
                }
            }

            // calculate average position
            float inv = 1.0f/(float)atomCounter;
            unsigned int keeper=0;
            xpos *= inv;
            ypos *= inv;
            zpos *= inv;
            // find in bead universe the bead that is closest
            float base = pModel->getBeadRadius()*1.13f;
            for(unsigned int b=0; b < totalBeads; b++){ // iterate over each bead in Universe
                currentBead = pModel->getBead(b);
                diffx = currentBead->getX() - xpos;
                diffy = currentBead->getY() - ypos;
                diffz = currentBead->getZ() - zpos;
                dis2 =std::sqrt((diffx*diffx + diffy*diffy + diffz*diffz));

                if (dis2 <= base){ //min = dis2;
                    std::cout << " => CENTERED BEAD FOUND " << b << " " << " RESID " << it->getResidByIndex(r_i) << std::endl;
                    keeper = b;
                    base = dis2;
                }
            }
            it->addCenteredAnchors(keeper);
        }
    }

    const float b2 = pModel->getBeadRadius()*1.13f;
    for(auto it = components.begin(); it != components.end(); ++it) {
        //Component temp = *it;

        for(unsigned int r=0; r<it->getTotalResids(); r++){
            std::cout << " SEARCHING ANCHOR " << it->getResidByIndex(r) << std::endl;
            for (unsigned int i=0; i < totalAtoms; i++){ // find all atoms that match resid and chain
                // match chain and resid to Component
                if ( it->getResidByIndex(r) == *(pdbResIDs + i) && (it->getChainsByIndex(r) == *(pdbChainIds + i)) ) {
                    xpos = *(pdbModel.getCenteredXVec() + i);
                    ypos = *(pdbModel.getCenteredYVec() + i);
                    zpos = *(pdbModel.getCenteredZVec() + i);
                    // find bead that is within radii
                    for(unsigned int b=0; b < totalBeads; b++){ // iterate over each bead in Universe
                        currentBead = pModel->getBead(b);
                        diffx = currentBead->getX() - xpos;
                        diffy = currentBead->getY() - ypos;
                        diffz = currentBead->getZ() - zpos;
                        float dis2 = std::sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

                        if (dis2 <= b2){
                            std::cout << " => ANCHOR ATOM FOUND " << pdbModel.getAtomTypeByIndex(i) << " " << *(pdbResIDs + i) << std::endl;
                            it->addAnchor(b);
                            break;
                        }
                    }
                }
            }
        }
    }

    anchorFile.close();
    anchorFile.clear();

    totalComponents = (unsigned int)components.size();

    std::cout << " TOTAL COMPONENTS " << totalComponents << std::endl;
//    for(auto it : components) {
//        auto sit =  it.getCenteredAnchors();
//        for (auto item : *sit){
//            std::cout << it.getID() << " " << item << std::endl;
//            Bead * pBead = pModel->getBead(item);
//            pBead->printCoordinates();
//        }
//    }

    if (components.empty())
        return false;

    return true;
}


/**
 * modelPR and targetPR are the same size
 * targetPR is derived from PDB
 */
float Anneal::calculateKLDivergenceAgainstPDBPR(std::vector<unsigned int> &modelPR, std::vector<double> &targetPR){

    double totalCounts = 0.0;
    double kl=0.0;
    double prob, * value;

    unsigned int totalm = modelPR.size();
    std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
    // normalization constant of model Pr
    // treats each value as discrete (i.e. not integrating via trapezoid)
    //last nonzero bin
    unsigned int last=0;
    for (unsigned int i=0; i < totalm; i++){
        if (targetPR[last] <= 0){
            break;
        }
        last++;
    }

    for (unsigned int i=0; i < totalm; i++){
        totalCounts += modelPR_float[i];
    }

    // for modelPR values in bins > shannon_bins are zero since p*log p/q = 0 for p=0
    for (unsigned int i=0; i < last; i++){
        prob = targetPR[i];  // get experimental bin value
        value = &modelPR_float[i]; // get model value
        if (prob > 0 && *value > 0){
            kl += prob * log(prob/(*value) * totalCounts);
        } else if (prob > 0 && *value <= 0){ // severely penalize any model bin that is zero
            kl += 10.0;
        }
    }

    return (float)(kl/(double)last);  // returns value per bin
}

/**
 * use to update distribution model
 * Input distribution must contain zereos for values not specified
 *
 * @param contactsDistributionSeed
 */
void Anneal::setContactsDistribution(std::vector<unsigned int> & contactsDistributionSeed){

    distributionlimit=contactsDistributionSeed.size();
    contactsDistribution.resize(distributionlimit);

    std::copy(contactsDistributionSeed.begin(), contactsDistributionSeed.end(), contactsDistribution.begin());

    double totaltemp=0;
    for(auto & item : contactsDistribution){
        totaltemp += item;
    }

    for(unsigned int i=0; i<distributionlimit; i++){ // normalize
        contactsDistribution[i] *= 1.0/totaltemp;
    }
}





double Anneal::estimateMagnitudeOfDifferenceContactsPotential(unsigned int workingLimit,
                                                              std::vector<unsigned int> &bead_indices,
                                                              std::vector<unsigned int> &binCount,
                                                              PointSetModel *pModel) {

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<unsigned int> backUpState(bead_indices);

    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    const unsigned int noNeigborIndex = pModel->getNeighborLimit();

    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class
    std::vector<unsigned int> binCountBackUp(binCount);
//    std::vector<unsigned int> binCountBackUp(binCount.size());  // smallish vector, typically < 50
//    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

    double deltaContacts = 0.0, currentContacts = calculateAvgContactsPotential(binCount, workingLimit);
    for(unsigned int i=0; i<101; i++){

        if (distribution(gen) < 0.5){ // remove
            unsigned int original = bead_indices[randomIndex(gen)];
            removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);

            deltaContacts += std::abs(currentContacts - calculateAvgContactsPotential(binCount, workingLimit));

            auto beginIt = bead_indices.begin();
            auto beginBinCount = binCount.begin();
            restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);

        } else { //add

            unsigned int addMe = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            while(addMe == noNeigborIndex){
                addMe = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            }

            auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), addMe);
            // make the swap at the border (workingLimit)
            addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);  // alters backUpState

            addToPr(addMe, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            deltaContacts += std::abs(currentContacts - calculateAvgContactsPotential(binCount, workingLimit));

            auto beginBinCount = binCount.begin();
            restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
            beads_in_use_tree.erase(addMe);
        }

    }

    return deltaContacts/101.0d;
}
