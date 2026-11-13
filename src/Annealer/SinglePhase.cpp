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

#include "../Anneal.h"
#include <sstream>
#include <stdexcept>
#include "../HullVolumeCache.h"

/*
 * Speed changes in this file (same moves, same random draws, same acceptance tests as before):
 *
 *  - ContactTracker replaces beads_in_use_tree (std::set) for occupancy, contact numbers and the contact energy:
 *    O(12) per move instead of a full recount with std::set look-ups.
 *  - SortedSelection keeps bead_indices exactly as before (sorted prefix) but finds free lattice points in O(1)
 *    and undoes rejected moves instead of copying the whole vector from backUpState.
 *  - createInitialModelCVXHull: the convex hull volume is recomputed only when the hull can change
 *    (an added point outside the hull, a removed point that is a hull vertex); otherwise the current volume is used.
 */



bool Anneal::createInitialModelCVXHull(PointSetModel *pModel, PofRData *pData, std::string name) {
    std::cout << " --------------------------------------------------- " << std::endl;
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    contactCutOff = interconnectivityCutOff;
    // convert distances to ShannonBin membership
    maxbin = pModel->populateBins(pData);
    dmModel = pModel; dmData = pData;     // for direct mode (bins not stored)

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");
    /*
     * pBin holds a vector of small ints that are used to determine which bins to add to
     */
    unsigned short int * pBin = pModel->getPointerToBins(); // initialized as empty in Model class

    maxbin += 1;  // maximum bin index, so vector size must be +1

    logger("INITIAL MAX BIN",formatNumber(maxbin));
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", formatNumber(totalBins));
    logger("MAX MODEL N_S BINS",formatNumber(maxbin));
    logger("BINWIDTH",formatNumber(pData->getBinWidth(),3));
    logger("BEAD RADIUS",formatNumber(pModel->getBeadRadius(),3));

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's

    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    pModel->writeSubModelToFile(0, totalBeadsInSphere, bead_indices, "universe");

    std::mt19937 gen(seedForRun(0));

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.39);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    float radius_larger = pModel->getBeadRadius()*std::sqrt(7.0f)/2.0f;
    auto lowerN = (unsigned int)std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
    auto upperN = (unsigned int)std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
    // the bead count comes from the Porod volume in the P(r) file; a missing or zero VOLUME gives no beads, and qhull
    // then fails on an empty point set (segmentation fault)
    if (pData->getVolume() <= 0 || upperN < 4 || upperN > totalBeadsInSphere){
        std::ostringstream msg;
        msg << "\n  ERROR createInitialModelCVXHull: cannot size the initial model from the P(r) file " << pData->getFilename()
                  << "\n        VOLUME : " << pData->getVolume() << " A^3 gives " << upperN << " beads of radius " << pModel->getBeadRadius()
                  << " (need 4 to " << totalBeadsInSphere << ").\n        Set a Porod volume (VOLUME) in the P(r) file.";
        throw std::invalid_argument(msg.str());
    }
    // pick random number between lowerN and upperN
    unsigned int workingLimit = upperN;//number_of_beads_to_use(gen);
    /*
     * create initial model
     * randomize bead indices
     * sort to workingLimit
     */
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    char flags[] = "qhull FA";
    HullVolumeCache hullVolume(pModel, flags, bead_indices, workingLimit);
    float test_volume, current_volume = hullVolume.currentVolume();
    std::set<unsigned int> hull;  // sort hull points into bead_indices

    recalculateDeadLimit(workingLimit, bead_indices, hull, pModel);
    std::string switched = "HULL";

    // occupancy (was beads_in_use_tree) and bookkeeping of bead_indices
    ContactTracker occupancy(pModel->getDirectPointerToNeighborhood(), pModel->getSizeOfNeighborhood(), pModel->getNeighborLimit(), totalBeadsInSphere);
    occupancy.reset(bead_indices.begin(), workingLimit);
    SortedSelection selection(bead_indices, workingLimit, totalBeadsInSphere);

    // make copy as initial best model
    unsigned int lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

    logger("WORKINGLIMIT SET", formatNumber(workingLimit));
    // setup parameters for hull
    std::vector<coordT> points(3*(upperN+(int)(0.30*upperN)) + 3*(size_t)totalBeadsInSphere);
    float lowest_volume = current_volume;
    float initialVolume = current_volume;

    logger("CVX VOLUME", formatNumber(current_volume, 1));

    EulerTour eulerTour(bead_indices, workingLimit, pModel);
    unsigned int tempNumberOfComponents, currentNumberOfComponents  = eulerTour.getNumberOfComponents();
    logger("EULER TOUR(s)", formatNumber(currentNumberOfComponents, 1));

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel);

    auto currentKL = (float)pData->getScore(binCount);
    float lowest_kl = currentKL;
    float hlambda = 0.0001f;
    double muPercent = 0.1;
    double muConstant = currentKL*muPercent/((1.0 - muPercent)*current_volume/(double)workingLimit);

    float temp_volume_energy, current_volume_energy = (float)(muConstant*current_volume/(double)workingLimit);

    float testKL, test_energy, current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + current_volume_energy;

    logger("INITIAL SCORE", formatNumber(currentKL, 7));
    float lowest_energy = current_energy;

    auto lowTempStop =  (double)0.001;//highT;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;
    unsigned int failures=0;

    double sumIt = 0.0, counter=0.0;
    bool modulo = false, updateCVX = false;

    logger("STARTING CONSTANT TEMP SEARCH", formatNumber(currentNumberOfComponents));

    unsigned int high, hullsize=hull.size();
    unsigned int testIndex, rI;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    double probabilityAddRemove = 0.361;

    highTempRounds = (unsigned int)(totalBeadsInSphere/2*std::log(lowerN)*2.7);
    if (highTempRounds < 41379){
        highTempRounds = 41379;
    }

    unsigned int updater = 397;
    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        if (distribution(gen) < probabilityAddRemove ){

            if ( (distribution(gen) < 0.53 && workingLimit < upperN) || workingLimit < lowerN ){
                // ADD
                testIndex = getUseableNeighborWeighted(occupancy, pModel, bead_indices[randomIndex(gen)], gen);
                while(testIndex == pModel->getNeighborLimit()){
                    testIndex = getUseableNeighborWeighted(occupancy, pModel, bead_indices[randomIndex(gen)], gen);
                }

                selection.add(testIndex);                    // == find in tail, iter_swap, wl++, sort
                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = (float)pData->getScore(binCount);

                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                test_volume = hullVolume.afterAdd(testIndex, bead_indices, workingLimit);
                temp_volume_energy = (float)(muConstant*test_volume/(double)workingLimit);

                test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);

                if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                    occupancy.add(testIndex);
                    hullVolume.accept();
                    updateCVX = true;
                } else { // undo changes (rejecting)
                    selection.undoAdd(testIndex);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    hullVolume.reject();
                    eulerTour.removeNode(testIndex);
                }

            } else {
                // REMOVE: grab from randomized active_indices list
                rI = randomIndex(gen);
                testIndex = bead_indices[rI];

                removeFromPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                selection.removeAt(rI);                      // == wl--, iter_swap, sort
                testKL = (float)pData->getScore(binCount);

                test_volume = hullVolume.afterRemove(testIndex, bead_indices, workingLimit);
                tempNumberOfComponents = eulerTour.removeNode(testIndex);
                temp_volume_energy = (float)(muConstant*test_volume/(double)workingLimit);

                test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);
                if (test_energy < current_energy || (std::exp((current_energy - test_energy)*inv_kb_temp) > distribution(gen))) {
                    occupancy.remove(testIndex);
                    hullVolume.accept();
                    updateCVX = true;
                } else { // undo changes and move to next bead (rejecting)
                    eulerTour.addNode(testIndex, pModel);
                    selection.undoRemove();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    hullVolume.reject();
                }
            }

        } else {
            // POSITIONAL
            unsigned int swap1 = bead_indices[randomIndex(gen)];

            if (currentNumberOfComponents > 1 && modulo){ // only move CVX points if true
                std::vector<unsigned int> active_indices(workingLimit);
                unsigned int * const pSelections = bead_indices.data(); // initialized as empty in Model class
                for (unsigned int i = 0; i < workingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                    active_indices[i] = pSelections[i];
                }

                // needs to be optimized
                qh_new_qhull(3, workingLimit, points.data(), 0, flags, nullptr, nullptr);
                vertexT * vertices = qh vertex_list;
                auto totalV = (unsigned int)qh num_vertices;

                // only move CVX hull points
                std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
                for (unsigned int v = 0; v < totalV; v++) { //
                    indices_to_check[v] = active_indices[qh_pointid( vertices->point)];
                    vertices = vertices->next;
                }

                qh_freeqhull(true);
                std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
                swap1 = indices_to_check[randomVertices(gen)];
                modulo=false;
            } else {
                // move a high dense bead
                while((getProbabilityOfSelectedIndex(occupancy, swap1)) > distribution(gen)){
                    swap1 = bead_indices[randomIndex(gen)];
                }
            }

            // find bead to swap in active set
            const unsigned int position = selection.find(swap1);
            // remove selected index from P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // find better position
            occupancy.remove(swap1);
            eulerTour.removeNode(swap1);

            unsigned int neighbor = getUseableNeighborWeighted(occupancy, pModel, bead_indices[randomIndex(gen)], gen);
            while(neighbor == pModel->getNeighborLimit() || neighbor == swap1){
                neighbor = getUseableNeighborWeighted(occupancy, pModel, bead_indices[randomIndex(gen)], gen);
            }

            selection.swapAt(position, neighbor);            // == iter_swap with tail position, sort

            addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = (float)pData->getScore(binCount);
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            test_volume = hullVolume.afterSwap(swap1, neighbor, bead_indices, workingLimit);
            temp_volume_energy = (float)(muConstant*test_volume/(double)workingLimit);

            test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);

            if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) { //
                occupancy.add(neighbor);
                hullVolume.accept();
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                selection.undoSwap(swap1, neighbor);
                hullVolume.reject();
                eulerTour.removeNode(neighbor);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                occupancy.add(swap1);
            }
        }

        // Adaptive simulated annealing part
        if (updateCVX){
            currentNumberOfComponents = tempNumberOfComponents;
            currentKL = testKL;
            current_volume = test_volume;
            current_volume_energy = temp_volume_energy;
            current_energy = test_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy

            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

            if (current_volume_energy > 1.5*currentKL && currentNumberOfComponents == 1){
                muConstant = currentKL*muPercent/((1.0 - muPercent)*current_volume/(double)workingLimit);
                current_volume_energy = (float)(muConstant*current_volume/(double)workingLimit);
            }

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASAConstantTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        if (high%updater == 0){
            std::cout << "*******************                                        *******************" << std::endl;
            printf("       ACCEPTRATE : %5.3f       STEP : %-7i (%7i)\n", acceptRate, high, highTempRounds);
            printf("            GRAPH : %5i   FAILURES : %i %s\n", currentNumberOfComponents, failures, switched.c_str());
            printf("            LIMIT : %5i   UPPER >= %i LOWER <= %i POOL : %i\n", workingLimit, upperN, lowerN, hullsize);
            printf("           VOLUME : %.0f (%d) MU : %.2E \n", current_volume, upperV, muConstant);
            printf("         %s : %.4E E_VOL : %.2E ENRGY : %.4E \n", score, currentKL, current_volume_energy, current_energy);
            std::cout << "*******************                                        *******************" << std::endl;
        }
        if (high % 171 == 0) modulo = true;

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            lowestWorkingLimit = workingLimit;
            lowest_volume = current_volume;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            lowest_energy = current_energy;
            lowest_kl = currentKL;
        }
    }

    logger("HULL VOLUMES COMPUTED", formatNumber(hullVolume.recomputed));
    logger("HULL VOLUMES REUSED", formatNumber(hullVolume.reused));

    highTempStartForCooling = (float)lowTempStop;

    EulerTour tempEulerTour(lowest_bead_indices.begin(), lowestWorkingLimit, pModel);
    std::cout << " Lowest Tour : " << tempEulerTour.getNumberOfComponents() << std::endl;

    std::set<unsigned int> centeringhullpts;
    std::set<unsigned int> lowestSet(lowest_bead_indices.begin(), lowest_bead_indices.begin()+lowestWorkingLimit);
    getHullPoints(centeringhullpts, lowestSet, pModel);

    pModel->centerLatticeModel(&lowestWorkingLimit, lowest_bead_indices, centeringhullpts);

    tempEulerTour = EulerTour(lowest_bead_indices.begin(), lowestWorkingLimit, pModel);
    std::cout << " Lowest Tour AFTER  : " << tempEulerTour.getNumberOfComponents() << std::endl;

    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);

    pModel->setBeadAverageAndStdev(lowestWorkingLimit, lowestWorkingLimit*0.2f);

    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    printf("          AVERAGE => %.0f (%.0f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    printf("      TRK AVERAGE => %.2f (%i) \n", (counter > 0 ? sumIt/counter : 0.0), (int)counter);
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                 VOLUME                 *******************" << std::endl;
    printf("          INITIAL => %.0f \n", initialVolume);
    printf("           LOWEST => %.0f \n", lowest_volume);
    printf("            FINAL => %.0f\n", current_volume);
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                    KL                  *******************" << std::endl;
    printf("           LOWEST => %.4E \n", lowest_kl);

    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "    SHANNON BINS IN DATA : " << pData->getShannonBins() << std::endl;

    if (pData->getIsPr()){
        std::cout << "            ZERO BIN AT : " << pData->getZeroBin() << std::endl;
    }

    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << std::endl;
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        tempAverageContacts += occupancy.contacts(bead_indices[i]);
    }
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << std::endl;

    if (tempEulerTour.getNumberOfComponents() == 1 ) {
        pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);
        return true;
    } else {
        name += "_failed";
        pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, -g (current value => " << highTempRounds << " )" << std::endl;
        return true;
    }
}


std::string Anneal::refineHomogenousBodyASAHybridEx(PointSetModel *pModel, PofRData *pData, std::string outputname) {
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;

    float runPercentAddRemove = percentAddRemove;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000001 : (double)highTempStartForCooling;
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    const char * pOutputFileName = outputname.c_str() ;
    std::string status = "RAMPING TO TEMP";

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    std::mt19937 gen(seedForRun(1));
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's

    // copy Starting_Set from initial model
    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit();
    const auto lowerLimit = (unsigned int) (workingLimit*0.67);

    // reset iterators to internal bead_indices
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    // set deadLimit of the selected set
    // convert distances in Search Sphere to ShannonBin membership
    this->fillPrBinsAndAssignTotalBin(pModel, pData); // convert all paired-distances in universe to P(r) bins
    unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as empty in Model class

    /*
     * create set of beads within CONVEX HULL
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    std::string hullName = "hull";
    pModel->writeSetToFile(hull, hullName);

    /*
     * create cross-validation set or working observed probability distribution that encompasses search space
     */
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", formatNumber(totalBins));
    logger("MAX MODEL N_S BINS",formatNumber(maxbin));
    logger("BINWIDTH",formatNumber(pData->getBinWidth(),3));
    logger("BEAD RADIUS",formatNumber(pModel->getBeadRadius(),3));

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel);
    double testKL, currentKL = pData->getScore(binCount);
    double lowestKL = currentKL, startingKL = currentKL;

    // CVX hull calculation
    char flags[] = "qhull FA";
    float starting_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // make backup of states
    std::vector<unsigned int> binCountBackUp(binCount);  // smallish vector, typically < 50

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    // coupon collector's problem
    unsigned int updateCount = ccmultiple*getCoupons(hull.size());
    float step_limit = updateCount;

    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();
    unsigned int tempNumberOfComponents = currentNumberOfComponents;

    std::cout << "             EULER TOURS : " << currentNumberOfComponents << std::endl;

    bool isUpdated = false, once = true;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    double runtime, inv_kb_temp = 1.0/lowTempStop;

    char addRemoveText[75];

    // occupancy, contact numbers and contact energy (was beads_in_use_tree + calculateNumberOfPairwiseContacts)
    ContactTracker contacts(pModel->getDirectPointerToNeighborhood(), pModel->getSizeOfNeighborhood(), noNeigborIndex, totalBeadsInSphere);
    contacts.reset(bead_indices.begin(), workingLimit);
    SortedSelection selection(bead_indices, workingLimit, totalBeadsInSphere);

    // contacts potential
    double tempContactEnergy;
    double contactsEnergy = contacts.energy();
    double startingContactEnergy = contactsEnergy;

    /*
     * there may be a chance that currentKLDivContacts comes back as 0, so if using to set a constant like dividing by,
     * then there will be a chance of return NAN
     */
    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);

    double contactDistributionWeight = 0.05*currentKL/contactsEnergy;
    // calculate and set starting energy of the model
    double current_energy = currentKL +
            lambda*connectivityPotential(currentNumberOfComponents) +
            contactDistributionWeight*contactsEnergy;

    double this_energy, lowest_energy = current_energy, starting_energy=current_energy;

    std::clock_t startTime;
    int attempts=0, failures=0, updated=0;

    std::vector<unsigned int> selections(totalBeadsInSphere);
    unsigned int * const pSelections = &selections[0]; // initialized as empty in Model class
    for(unsigned int i=0; i < workingLimit; i++){
        *(pSelections+i) = i;
    }

    // before the loop
    WeightSchedule contactShare{{ {0.00, 0.05}, {0.15, 0.30}, {0.65, 0.30}, {1.00, 0.50} }};  // example
    double dKLmean = 0.0, dCmean = 0.0;
    const unsigned int warmup = 200;                                 // proposals before the first weight update
    const unsigned int updateEvery = std::max(1u, (unsigned int)(0.005*step_limit));   // every 0.5% of the run

    unsigned int compCOunt=0, numberOfCoolingTempSteps=0;
    unsigned int addremoveCount = 0, moveCount = 0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        startTime = std::clock();

        if ( distribution(gen) < runPercentAddRemove) { //add or remove bead within working Set (exclude deadzone)
            addremoveCount++;
            if (workingLimit < lowerLimit || distribution(gen) < 0.493 ){ // ADD BEAD?
                /*
                 * pick random bead to ADD TO.  treats each bead equally which means that the neighbor we pick is biased
                 * to denser regions - only look at available neighbors so the new model is fully connected
                 */
                unsigned int addMe = getUseableNeighborFromSet(contacts, pModel, bead_indices[randomIndex(gen)], gen);
                while(addMe == noNeigborIndex){
                    addMe = getUseableNeighborFromSet(contacts, pModel, bead_indices[randomIndex(gen)], gen);
                }

                selection.add(addMe);                        // == find in tail, iter_swap, wl++, sort
                addToPr(addMe, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->getScore(binCount);
                contacts.add(addMe);
                tempContactEnergy = contacts.energy();

                tempNumberOfComponents = eulerTour.addNode(addMe, pModel);

                this_energy = testKL +
                    lambda*connectivityPotential(tempNumberOfComponents) +
                        contactDistributionWeight*tempContactEnergy;

                if ( (this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) ))) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "   ADDED => %i", addMe);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1);
                } else { // undo changes (rejecting)
                    contacts.remove(addMe);
                    eulerTour.removeNode(addMe);
                    std::sprintf(addRemoveText, "     ADD => %i FAILED", addMe);
                    selection.undoAdd(addMe);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                }

            } else { // REMOVE BEADS?
                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                unsigned int original = bead_indices[selections[0]];

                bool foundIt = true;
                for(unsigned int s=1; s<workingLimit; s++){ // randomly remove bead constrained by number of components
                    tempNumberOfComponents = eulerTour.removeNode(original);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        foundIt = false;
                        break;
                    }
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[selections[s]];
                }

                if (foundIt){
                    original = bead_indices[randomIndex(gen)];
                    tempNumberOfComponents = eulerTour.removeNode(original);
                }

                const unsigned int k = selection.find(original);
                removeFromPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                selection.removeAt(k);                       // == wl--, iter_swap, sort

                testKL = pData->getScore(binCount);

                contacts.remove(original);
                tempContactEnergy = contacts.energy();

                this_energy = testKL +
                        lambda*connectivityPotential(tempNumberOfComponents) +
                        contactDistributionWeight*tempContactEnergy;

                if ( ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) ))) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, " REMOVED => %i", original);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1);
                } else { // undo changes and move to next bead (rejecting)
                    contacts.add(original);
                    selection.undoRemove();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());

                    eulerTour.addNode(original, pModel);
                    std::sprintf(addRemoveText, "  REMOVE => %i FAILED", original);
                }
            }

        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end
            attempts +=1;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            moveCount++;
            // pick a point to move, we want to favor denser regions (i.e., points with high contact numbers)
            bool foundIt = true;

            // only remove bead that is singly connected
            if (numberOfCoolingTempSteps > sixtyfivePercent) {
                for(unsigned int s=0; s<workingLimit; s++){
                    position = selections[s];
                    swap1 = bead_indices[position];
                    if (contacts.contacts(swap1) < 2 ) {
                        // favors beads with greater contacts
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents) {
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if (foundIt && numberOfCoolingTempSteps > sixtyfivePercent) {
                for(unsigned int s=0; s<workingLimit; s++){
                    position = selections[s];
                    swap1 = bead_indices[position];
                    auto tt = contacts.contacts(swap1);
                    if (tt < 5 || tt > 8) {
                        // favors beads with greater contacts
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents) {
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if(foundIt){ // if nothing found because of probabilistic selection - then grab 1st best one
                for(unsigned int s=0; s<workingLimit; s++){
                    position = selections[s];
                    swap1 = bead_indices[position];
                    // favors beads with greater contacts
                    if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                        compCOunt++;
                        break;
                    }
                    eulerTour.addNode(swap1, pModel);
                }
            }

            // remove selected index from P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

            contacts.remove(swap1);

            // find new position weighted by connectivity, fewer connections has higher probability
            unsigned int originalSwap2Value = noNeigborIndex; // find new position weighted by connectivity
            for(unsigned int s=0; s<workingLimit; s++){ // array indices max value is workingLimit
                unsigned int selected = bead_indices[selections[s]];
                originalSwap2Value = getUseableNeighborFromSet(contacts, pModel, selected, gen);
                if (originalSwap2Value != noNeigborIndex && originalSwap2Value != swap1 && selected != swap1){
                    break;
                }
            }

            if (originalSwap2Value == noNeigborIndex || originalSwap2Value == swap1){
                // no free position found (previously undefined behaviour: iter_swap with end()); reject the move
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                contacts.add(swap1);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                std::sprintf(addRemoveText, "  FAILED => %i", swap1);
                tempContactEnergy = contactsEnergy;
                testKL = currentKL;
            } else {

                // make the swap, sort and update P(r)
                selection.swapAt(position, originalSwap2Value);  // == iter_swap with tail position, sort

                addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                testKL = pData->getScore(binCount);

                contacts.add(originalSwap2Value); // add new lattice
                tempContactEnergy = contacts.energy();

                this_energy = testKL + contactDistributionWeight*tempContactEnergy;
                tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
                this_energy += lambda*connectivityPotential(tempNumberOfComponents);

                if ( (this_energy < current_energy || (std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)))) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, " SWAPPED => %i to %i ", swap1, originalSwap2Value);
                } else {
                    selection.undoSwap(swap1, originalSwap2Value);
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                    contacts.remove(originalSwap2Value);
                    contacts.add(swap1); // add back swap one
                    std::sprintf(addRemoveText, "  FAILED => %i", swap1);
                    eulerTour.removeNode(originalSwap2Value);
                    currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                }
            }
        } // end of positional refinement or add/remove if statement

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        // typical size of each term's change per proposal (accepted or not)
        dKLmean = 0.99*dKLmean + 0.01*std::fabs(testKL - currentKL);
        dCmean  = 0.99*dCmean  + 0.01*std::fabs(tempContactEnergy - contactsEnergy);

        // Adaptive simulated annealing part
        if (isUpdated){
            currentKL = testKL;
            currentNumberOfComponents = tempNumberOfComponents;
            contactsEnergy = tempContactEnergy;

            current_energy = this_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i;
            }

            if (numberOfCoolingTempSteps > sixtyfivePercent){
                status = once ? "ANNEALING" : status;

                if (current_energy < lowest_energy){
                    lowestKL = currentKL;
                    lowest_energy = current_energy;
                }
            }

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        // phase labels only (no weight changes here)
        if (numberOfCoolingTempSteps < fifteenPercent)        status = "RAMPING";
        else if (numberOfCoolingTempSteps < sixtyfivePercent) status = "EQUILIBRATING";

        // contact weight from the schedule: contacts carry share s of the typical energy change per move
        if (numberOfCoolingTempSteps >= warmup && numberOfCoolingTempSteps % updateEvery == 0 && dCmean > 0 && dKLmean > 0) {
            const double s = contactShare.at((double)numberOfCoolingTempSteps / step_limit);
            contactDistributionWeight = s / (1.0 - s) * dKLmean / dCmean;
            current_energy = currentKL + contactDistributionWeight*contactsEnergy
                           + lambda*connectivityPotential(currentNumberOfComponents);
            lowest_energy = current_energy;          // energies under different weights are not comparable
        }

        if (numberOfCoolingTempSteps == sixtyfivePercent){
            std::string name = "equilibrated";
            pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, name, this, pData,
                                      numberOfCoolingTempSteps, 120000, contactsEnergy);
        }

        if (numberOfCoolingTempSteps % 1001 == 0) {
            float percentADDREMOVEPos = (float)addremoveCount/(float)(addremoveCount+moveCount);

            printf("      TEMP : %-.3E MAXSTEPS => %.0f (%i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, status.c_str());
            printf("    ACCEPT : %9.5f FAILURES => %i  %.2f TIME : %.4f\n", acceptRate, failures, percentADDREMOVEPos, runtime);
            printf("     LIMIT : %9i  (LL %6i) COMP => %3i (%i) \n", workingLimit, lowerLimit, currentNumberOfComponents, compCOunt);
            printf("  CONTACTS : %.3E (%.1E | %.1E) \n", contactDistributionWeight*tempContactEnergy, contactsEnergy, contactDistributionWeight);
            printf("           : %s %s\n", pOutputFileName, addRemoveText);
            printf("   %s => %-4.3E ( %.4E ) ENRGY : %.4E ( %.4E )\n", score, currentKL, lowestKL, current_energy, starting_energy);
            moveCount=0;
            addremoveCount=0;
        }

        // random access
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = currentKL;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = contactsEnergy;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

#ifdef ANNEAL_VERIFY
        // incremental state against full recomputation (debug builds only)
        {
            std::vector<unsigned int> fullBins(maxbin);
            calculateModelPrDistribution(&bead_indices, &fullBins, workingLimit, totalBeadsInSphere, pModel);
            std::set<unsigned int> tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
            if (fullBins != binCount || calculateNumberOfPairwiseContacts(workingLimit, pModel, tree, bead_indices) != contactsEnergy
                || !std::is_sorted(bead_indices.begin(), bead_indices.begin() + workingLimit) || eulerTour.validateList("refine")){
                std::cout << "VERIFY FAILED at step " << numberOfCoolingTempSteps << std::endl;
                exit(1);
            }
        }
#endif

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
    } // end of steps

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    double totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (contacts.contacts(bead_indices[i]) == c){
                totalCounts += 1.0;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (unsigned int c=1; c<13; c++){
        double totalContactsAt = 0.0;
        for (unsigned int i=0; i<workingLimit; i++){
            if (contacts.contacts(bead_indices[i]) == c){
                totalContactsAt += 1.0;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/totalCounts);
    }

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    pModel->setCVXHullVolume(current_volume);
    // reset temp for low temp annealing (strictly positional)
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        tempAverageContacts += contacts.contacts(bead_indices[i]);
    }

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    std::cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << std::endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << " FINAL LOW TEMP REFINEMENT" << std::endl;
    std::cout << "              KL DIVERGENCE   ENERGY " << std::endl;
    std::printf("      START =>    %.4E           \n", startingKL);
    std::printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    std::cout << "    CONTACTS" << std::endl;
    std::printf("      AVERAGE =>    %.2f TOTAL => %.0f \n", average_number_of_contacts, tempAverageContacts);
    std::printf("    START => %.4E \n", startingContactEnergy);
    std::printf("    FINAL => %.4E \n", contactsEnergy);
    std::printf("    VOLUME  START => %.4E FINAL => %.4E \n", starting_volume, current_volume);
    std::cout << "------------------------------------------------------------------------------" << std::endl;

    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            outputname,
            this,
            pData,
            numberOfCoolingTempSteps,
            current_volume,
            average_number_of_contacts);

    return nameOfModel;
}
