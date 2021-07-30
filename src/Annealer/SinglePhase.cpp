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


bool Anneal::createInitialModelCVXHull(PointSetModel *pModel, PofRData *pData, std::string name) {
    std::cout << " --------------------------------------------------- " << std::endl;
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    contactCutOff = interconnectivityCutOff;
    // convert distances to ShannonBin membership
    maxbin = pModel->populateBins(pData);

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
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", formatNumber(totalBins));
    logger("MAX MODEL N_S BINS",formatNumber(maxbin));
    logger("BINWIDTH",formatNumber(pData->getBinWidth(),3));
    logger("BEAD RADIUS",formatNumber(pModel->getBeadRadius(),3));

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    // c-style is slightly faster for large vector sizes
    //std::clock_t start = std::clock();
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }
    // cout << "C-STYLE " << (std::clock() - start)/(double) CLOCKS_PER_SEC << endl;

    pModel->writeSubModelToFile(0, totalBeadsInSphere, bead_indices, "universe");


    std::random_device rd;
    std::mt19937 gen(rd());

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.39);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    float radius_larger = pModel->getBeadRadius()*std::sqrt(7.0f)/2.0f;
    auto lowerN = (unsigned int)std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
    auto upperN = (unsigned int)std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
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
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    std::set<unsigned int> hull;  // sort hull points into bead_indices
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    //populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    recalculateDeadLimit(workingLimit, bead_indices, hull, pModel);
    std::string switched = "HULL";

    // make copy as initial best model
    unsigned int lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

    logger("WORKINGLIMIT SET", formatNumber(workingLimit));
    // setup parameters for hull
    coordT points[3*(upperN+(int)(0.30*upperN))];
    //char flags[] = "qhull FA";
    //float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
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
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
//    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize); // guaranteed unbiased

    double probabilityAddRemove = 0.361;

    highTempRounds = (unsigned int)(totalBeadsInSphere/2*std::log(lowerN)*2.7);
    if (highTempRounds < 41379){
        highTempRounds = 41379;
    }

    unsigned int updater = 397;
    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        if (distribution(gen) < probabilityAddRemove ){

            if ( (distribution(gen) < 0.53 && workingLimit < upperN) || workingLimit < lowerN ){
//                std::cout << "*******************                  ADD                   ******************* "  << std::endl;

                testIndex = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while(testIndex == pModel->getNeighborLimit()){
                    testIndex = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }

//                auto set_it=hull.begin();
//                std::advance(set_it, randomHull(gen));
//                testIndex = *set_it;
//
//                while ( beads_in_use_tree.find(testIndex) != beads_in_use_tree.end()) { // find a new position
//                    //testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
//                    set_it=hull.begin();
//                    std::advance(set_it, randomHull(gen));
//                    testIndex = *set_it;
//                }

                auto itIndex = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), testIndex);
                //make the swap at the border (workingLimit)

                addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);
                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = (float)pData->getScore(binCount);

                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                // will be very slow calculation for large systems!
                test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                temp_volume_energy = (float)(muConstant*test_volume/(double)workingLimit);

                test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);

                if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                    beads_in_use_tree.insert(testIndex);
                    updateCVX = true;
                } else { // undo changes (rejecting)
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                }

            } else {
//                std::cout << "*******************                 REMOVE                 ******************* "  << std::endl;
                // grab from randomized active_indices list
                rI = randomIndex(gen);
                testIndex = bead_indices[rI];

                removeLatticePositionToModelByIndex(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, rI);
                testKL = (float)pData->getScore(binCount);

                test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                tempNumberOfComponents = eulerTour.removeNode(testIndex);
                temp_volume_energy = (float)(muConstant*test_volume/(double)workingLimit);

                test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);
                if (test_energy < current_energy || (std::exp((current_energy - test_energy)*inv_kb_temp) > distribution(gen))) {
                    beads_in_use_tree.erase(testIndex);
                    updateCVX = true;
                } else { // undo changes and move to next bead (rejecting)
                    eulerTour.addNode(testIndex, pModel);
                    auto beginIt = bead_indices.begin();
                    auto beginBinCount = binCount.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                }
            }


        } else {

//            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            // recalculate
            // calculate convex hull and get hull points
            // can be threaded
            unsigned int swap1 = bead_indices[randomIndex(gen)];

            if (currentNumberOfComponents > 1 && modulo){ // only move CVX points if true
                std::vector<unsigned int> active_indices(workingLimit);
                unsigned int * const pSelections = bead_indices.data(); // initialized as empty in Model class
                for (unsigned int i = 0; i < workingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                    active_indices[i] = pSelections[i];
                }

                // needs to be optimized
                qh_new_qhull(3, workingLimit, points, 0, flags, nullptr, nullptr);
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
                while((getProbabilityOfSelectedIndex(&beads_in_use_tree, pModel, swap1)) > distribution(gen)){
                    swap1 = bead_indices[randomIndex(gen)];
                }
            }

            // find bead to swap in active set
            auto itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
            // remove selected index from P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // find better position
            beads_in_use_tree.erase(swap1);
            eulerTour.removeNode(swap1);
            // find a bead in use, and see if it has an empty neighbor to use
            //unsigned int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]); // biased towards denser regions
            //while (neighbor == noNeigborIndex || neighbor == swap1){
            //    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            //}
            // if neighbor is not in use, then find it in bead_indices and assign to pSwap2
            // make the swap, sort and update P(r)

//            auto set_it=hull.begin();
//            std::advance(set_it, randomHull(gen));
//            unsigned int neighbor = *set_it;
//
//            while (beads_in_use_tree.find(neighbor) != beads_in_use_tree.end() || neighbor == swap1) { // find a new position
//                //testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
//                set_it=hull.begin();
//                std::advance(set_it, randomHull(gen));
//                neighbor = *set_it;
//            }


            unsigned int neighbor = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            while(neighbor == pModel->getNeighborLimit() || neighbor == swap1){
                neighbor = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            }

            auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), neighbor); // potential problem

            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = (float)pData->getScore(binCount);
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
            temp_volume_energy = (float)(muConstant*test_volume/(double)workingLimit);

            test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);

            if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) { //
                beads_in_use_tree.insert(neighbor);
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                eulerTour.removeNode(neighbor);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                beads_in_use_tree.insert(swap1);
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

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy

//            if (currentNumberOfComponents < 13){
//                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel); // list of points that are connected to neighbors
//                switched = "NEAREST NEIGHBOR";
//            } else {
//                switched = "HULL";
//                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel); // everything in CVX hull
//            }
//
//            hullsize = hull.size();
//            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);
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

        modulo = (high%171) ? true : modulo;
        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }
//
        // test bead_indices and vector are in sync
//        for(int i=0; i<workingLimit; i++){
//            if (beads_in_use_tree.find(bead_indices[i]) == beads_in_use_tree.end()){
//                cout << "!!!!! [ " << i << " ] NOT FOUND IN TREE => " << bead_indices[i] << endl;
//                exit(0);
//            }
//        }
//
//        for(auto bit : beads_in_use_tree){
//            if(std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, bit) == bead_indices.begin()+workingLimit){
//                std::cout << "!!!!! " << " NOT FOUND in BEAD_INDICES " << bit << std::endl;
//                exit(0);
//            }
//        }
//
//        if(beads_in_use_tree.size() != workingLimit){
//            cout << "!!!!! " << " Incorrect size WORKINGLIMIT "  << endl;
//            exit(0);
//        }

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            lowestWorkingLimit = workingLimit;
            lowest_volume = current_volume;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            lowest_energy = current_energy;
            lowest_kl = currentKL;
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }
    }

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
    printf("      TRK AVERAGE => %.2f (%i) \n", (sumIt/counter), (int)counter);
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

//    pData->writeICalc(binCount);
    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << std::endl;
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
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

    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000001 : (double)highTempStartForCooling;
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    const char * pOutputFileName = outputname.c_str() ;
    std::string status = "RAMPING TO TEMP";

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere); // large vector ~1000's
//    std::vector<unsigned int> backUpState(totalBeadsInSphere);

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
    this->fillPrBinsAndAssignTotalBin(pModel, pData);
    unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    /*
     * create set of beads within CONVEX HULL
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    std::string hullName = "hull";
    pModel->writeSetToFile(hull, hullName);
    std::set<unsigned int> original_hull(beads_in_use_tree.begin(), beads_in_use_tree.end());

    const auto targetCount = (unsigned int)(0.8*(original_hull.size()));
    /*
     * create cross-validation set or working observed probability distribution that encompasses search space
     */
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
//    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

    logger("TOTAL EXP N_S BINS", formatNumber(totalBins));
    logger("MAX MODEL N_S BINS",formatNumber(maxbin));
    logger("BINWIDTH",formatNumber(pData->getBinWidth(),3));
    logger("BEAD RADIUS",formatNumber(pModel->getBeadRadius(),3));

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel);
    double testKL, currentKL = pData->getScore(binCount);
    double lowestKL = currentKL, startingKL = currentKL;

    unsigned int exponent = std::ceil(-1.0*std::log10(currentKL));
    startingKL = std::ceil(std::pow(10,exponent)*currentKL)/std::pow(10,exponent);
    // CVX hull calculation
    char flags[] = "qhull FA";
    float starting_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    std::vector<unsigned int> binCountBackUp(binCount);  // smallish vector, typically < 50
//    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    /*
     * model for initial search is too dense, so deadlimit layer should be sufficient for a confined search
     */
    std::vector<unsigned int> backUpState(bead_indices);
    //std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    // coupon collector's problem
    unsigned int updateCount = ccmultiple*getCoupons(hull.size());

    // float step_limit = updateCount/0.65f;///0.1;
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

    bool isUpdated = false, once = true, writeOnce = true;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    double runtime, inv_kb_temp = 1.0/lowTempStop;
    //int diffContacts = pModel->getSizeOfNeighborhood() - this->contactsPerBead;
    char addRemoveText[75];

    // contacts potential
//    double tempContactSum, currentContactsSum = contactPotential(&beads_in_use_tree, pModel);
//    double tempKLDivContacts, currentKLDivContacts = currentContactsSum*invTotal;
    //float target = 9.0f;
    targetContacts = 1.9111;
//      targetContacts = 2.51111;
//    targetContacts = 2.7111;
//    targetContacts = 8.7111;
//    targetContacts = 4.11;
//    targetContacts = 0.99*contactsAvg;
    double contactsAvg = 2.0*binCount[0]/(float)workingLimit;
//    double currentContactsPotential = calculateAvgContactsPotential(binCount, workingLimit);

    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
    double tempKLDivContacts, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);
    /*
     * there may be a chance that currentKLDivContacts comes back as 0, so if using to set a constant like dividing by,
     * then there will be a chance of return NAN
     */

//    double deltaContacts = estimateMagnitudeOfDifferenceContactsPotential(workingLimit, bead_indices, binCount, pModel);
//    deltaContacts = (deltaContacts == 0) ? currentContactsPotential : deltaContacts;

    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);
    const auto eightyfivePercent = (unsigned int)(0.85*step_limit);
    const auto ninetyfivePercent = (unsigned int)(step_limit - lowerLimit/(1.0-percentAddRemove));

    // slowly increase weight of total Contact energy over D_kl
//    double currentCDW = alpha/1000;
//    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
//    double contactDistributionWeight = alphaConstant*currentKL/currentKLDivContacts;
//    double scaleFactor = alpha*currentKL/currentKLDivContacts;
    double targetKL = 0.00001;// should be relative to 10^-6
    double scaleFactor = alpha*targetKL/0.01;

//    double scaleFactor = alpha*targetKL/currentKLDivContacts;
//    double startingScaleReduction = 1000;
    double currentCDW = scaleFactor;///startingScaleReduction;
//    double alphaConstant = alpha > 0 ? std::pow(startingScaleReduction, 1.0/(double)fifteenPercent) : 0;
//    alphaConstant = 0;
    double contactDistributionWeight = currentCDW;
    //double contactDistributionWeight = currentKL*alphaConstant/((1.0 - currentCDW)*currentKLDivContacts);

    double current_energy = currentKL +
            lambda*connectivityPotential(currentNumberOfComponents) +
            contactDistributionWeight*currentKLDivContacts;

    double this_energy, lowest_energy = current_energy, starting_energy=current_energy;

    std::clock_t startTime;
    int attempts=0, failures=0, updated=0;

    std::vector<unsigned int> selections(totalBeadsInSphere);
    unsigned int * const pSelections = &selections[0]; // initialized as empty in Model class
    for(unsigned int i=0; i < workingLimit; i++){
        *(pSelections+i) = i;
    }

    unsigned int compCOunt=0, numberOfCoolingTempSteps=0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        startTime = std::clock();

        if ( distribution(gen) < percentAddRemove) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // build a list of indices within defined region of convex hull
            if (workingLimit < lowerLimit || distribution(gen) < 0.51 ){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                /*
                 * pick random bead that probablistically favors the least connected
                 */
//                unsigned int addMe = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
//                while(addMe == noNeigborIndex){
//                    addMe = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
//                }

                /*
                 * pick random bead to add to treats each bead equally which means that the neighbor we pick is biased
                 * to denser regions
                 */
                unsigned int addMe = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while(addMe == noNeigborIndex){
                    addMe = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }

                auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), addMe);

                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);  // alters backUpState
                addToPr(addMe, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->getScore(binCount);

                beads_in_use_tree.insert(addMe);
                addToContactsDistribution(addMe, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);

                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

                tempNumberOfComponents = eulerTour.addNode(addMe, pModel);

                this_energy = testKL +
                        lambda*connectivityPotential(tempNumberOfComponents) +
                        contactDistributionWeight*tempKLDivContacts;

                //testKL <= startingKL &&
                //this_energy <= starting_energy &&
                if ( (this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) ))) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "   ADDED => %i", addMe);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(addMe);
                    eulerTour.removeNode(addMe);
                    std::sprintf(addRemoveText, "     ADD => %i FAILED", addMe);
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
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

                removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);

                // still need to sort, swap changes the order
                testKL = pData->getScore(binCount);

                removeFromContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                beads_in_use_tree.erase(original);

                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
                this_energy = testKL +
                        lambda*connectivityPotential(tempNumberOfComponents) +
                        contactDistributionWeight*tempKLDivContacts;

                if ( ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) ))) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, " REMOVED => %i", original);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1);
                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
                    auto beginIt = bead_indices.begin();
                    auto beginBinCount = binCount.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);

                    eulerTour.addNode(original, pModel);
                    std::sprintf(addRemoveText, "  REMOVE => %i FAILED", original);
                }
            }

        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end
            attempts +=1;
            // shuffling should not change the location of the iterator
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            // pick a point to move, we want to favor denser regions (i.e., points with high contact numbers)
            bool retry = true;

            if (numberOfCoolingTempSteps > ninetyfivePercent){ // strictly favors points that are singly connected
                for(unsigned int s=1; s<workingLimit; s++){
                    swap1 = bead_indices[position];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1) < 2){
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            retry = false;
                            compCOunt++;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                    position = selections[s];
                }
            }

            if (retry){ // favors beads that are well connected
                bool foundIt = true;
//                for(unsigned int s=0; s<workingLimit; s++){
//                    position = selections[s];
//                    swap1 = bead_indices[position];
//                    if ((getProbabilityOfSelectedIndex(&beads_in_use_tree, pModel, swap1)) < distribution(gen)){ // favors beads with greater contacts
//                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
//                            compCOunt++;
//                            foundIt = false;
//                            break;
//                        }
//                        eulerTour.addNode(swap1, pModel);
//                    }
//                }

                if(foundIt){ // if nothing found because of probablistic selection - then grap 1st best one
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
            }

            // remove contribution of swap1
            auto itIndex = bead_indices.begin() + position;

            // remove selected index from P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

            // remove from beads_in_use_tree
            removeFromContactsDistribution(swap1, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            beads_in_use_tree.erase(swap1);

            // find new position weighted by connectivity, fewer connections has higher probability
            unsigned int originalSwap2Value = noNeigborIndex; // find new position weighted by connectivity
            for(unsigned int s=0; s<workingLimit; s++){ // array indices max value is workingLimit
                unsigned int selected = bead_indices[selections[s]];
                originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, selected);
                //originalSwap2Value = getUseableNeighborWeighted(&beads_in_use_tree, pModel, selected);
                if (originalSwap2Value != noNeigborIndex && originalSwap2Value != swap1 && selected != swap1){
                    break;
                }
            }

            // make the swap, sort and update P(r)
            auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);

            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            addToContactsDistribution(originalSwap2Value, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);

            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

            this_energy = testKL + contactDistributionWeight*tempKLDivContacts;
            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
            this_energy += lambda*connectivityPotential(tempNumberOfComponents);

            //testKL <= startingKL &&
            if ( (this_energy < current_energy || (std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)))) {
                isUpdated = true;
                std::sprintf(addRemoveText, " SWAPPED => %i to %i ", swap1, originalSwap2Value);
            } else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                beads_in_use_tree.erase(originalSwap2Value);
                std::sprintf(addRemoveText, "  FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            currentKLDivContacts = tempKLDivContacts;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

            //unsigned int * const pSelections = &selections[0]; // initialized as empty in Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i;
            }

            if (numberOfCoolingTempSteps > sixtyfivePercent){
                status = once ? "ANNEALING" : status;
//                if (once && numberOfCoolingTempSteps > ninetyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
//                    contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
//                    current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
//                    lowest_energy = current_energy;
//                    once = false;
//                    status = "FINAL LOW TEMP EQ";
//                }

                if (current_energy < lowest_energy){
                    lowestKL = currentKL;
                    lowest_energy = current_energy;
                }
            }

            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/(float)workingLimit;
        } else {
            std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        // rescale contactsWeight in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent){
//            currentCDW *= alphaConstant;
//            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
//            contactDistributionWeight *= alphaConstant;

            current_energy = (float)(currentKL +
                    lambda*connectivityPotential(currentNumberOfComponents) +
                    contactDistributionWeight*currentKLDivContacts);

            lowest_energy = current_energy;
            status = "RAMPING";
        } else if (numberOfCoolingTempSteps < sixtyfivePercent) {
            status = "EQUILIBRATING";
        }

        if (numberOfCoolingTempSteps > sixtyfivePercent && writeOnce){
                std::string name = "equilibrated";
                std::string nameOfModel = pModel->writeModelToFile2(
                        currentKL,
                        workingLimit,
                        bead_indices,
                        binCount,
                        name,
                        this,
                        pData,
                        numberOfCoolingTempSteps,
                        120000,
                        contactsAvg);
                writeOnce = false;
        }


        double tempDKLValue = contactDistributionWeight*currentKLDivContacts;

        printf("      TEMP : %-.3E MAXSTEPS => %.0f (%i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, status.c_str());
        printf("    ACCEPT : %9.5f FAILURES => %i  TIME : %.4f\n", acceptRate, failures, runtime);
        printf("     LIMIT : %9i  COMP => %3i (%i) \n", workingLimit, currentNumberOfComponents, compCOunt);
        printf("  CONTACTS : %.3E (%.1E | %.1E) ALPHA : %.1E AVG %.2f\n", tempDKLValue, currentKLDivContacts, contactDistributionWeight, contactDistributionWeight, contactsAvg);
        printf("           : %s %s\n", pOutputFileName, addRemoveText);
        printf("   %s => %-4.3E ( %.4E ) ENRGY : %.4E ( %.4E )\n", score, currentKL, lowestKL, current_energy, starting_energy);

        // random access
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = currentKL;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKLDivContacts;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        //calculateModelPrDistributionDirect(&bead_indices, &testBinCount, workingLimit, pModel, pData);
        //calculateModelPrDistribution(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        float testKL1 = pData->getScore(testBinCount);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            std::cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << std::endl;
//            return "stopped";
//        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
    } // end of steps

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    unsigned int originalCount=0;
    for(auto sit = beads_in_use_tree.begin(); sit != beads_in_use_tree.end(); ++sit){
        if (original_hull.find(*sit) != original_hull.end()){
            originalCount++;
        }
    }

    std::cout << " ORIGINAL FROM HULL TOTAL : " << originalCount << " target count " << targetCount << std::endl;
//    std::string name = "ar_" + outputname;
//    pModel->writeModelToFile2(
//            currentKL,
//            workingLimit,
//            bead_indices,
//            binCount,
//            name,
//            this,
//            pData,
//            numberOfCoolingTempSteps,
//            10000,
//            0);
//    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
//    currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    double totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts += 1.0d;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (unsigned int c=1; c<13; c++){
        double totalContactsAt = 0.0d;
        for (unsigned int i=0; i<workingLimit; i++){
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt += 1.0d;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/totalCounts);
    }

//    std::cout << "KL DIVERGENCE " << std::endl;
//    pData->printKLDivergence(binCount);

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    pModel->setCVXHullVolume(current_volume);
    // reset temp for low temp annealing (strictly positional)
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    // CONSTANT TEMP REFINEMENT
    // make Contact Potential a percentage of final energy
//    // At end of each temp, update a probability model for volume?  Use this to select
//    // string tempName = "rough";
//    // pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, tempName, this, pData);
//
//    pModel->writeModelToFile(workingLimit, bead_indices, "before_final", numberOfCoolingTempSteps);
//    cout << "------------------------------------------------------------------------------" << endl;
//    printf(" NUMBER OF STEPS %i\n", numberOfCoolingTempSteps);
//    printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
//    cout << "------------------------------------------------------------------------------" << endl;
//

    std::cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << std::endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << " FINAL LOW TEMP REFINEMENT" << std::endl;
    std::cout << "              KL DIVERGENCE   ENERGY " << std::endl;
    std::printf("      START =>    %.4E           \n", startingKL);
    std::printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    std::cout << "    CONTACTS" << std::endl;
    std::printf("      AVERAGE =>    %.2f TOTAL => %.0f \n", average_number_of_contacts, tempAverageContacts);
    std::printf("    BIN [ONE] => %d AVG => %.3f \n", binCount[0], binCount[0]/(float)workingLimit);
    std::printf("    VOLUME  START => %.4E FINAL => %.4E \n", starting_volume, current_volume);
    std::cout << "------------------------------------------------------------------------------" << std::endl;

//
//   this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);
//   //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
//
//    // if multithread, must put a lock on these two step
//    //CVX HULL STUFF
    /*
     * print distances
     */
//    for(unsigned int i=0; i<workingLimit; i++){
//
//        unsigned next = i+1;
//        const vector3 * vec1 = &pModel->getBead(bead_indices[i])->getVec();
//
//        for(; next<workingLimit; next++){
//            auto vec3 = *vec1 - pModel->getBead(bead_indices[next])->getVec();
//            std::cout << i << " " << next << " " << vec3.length() << " " << pData->convertToBin(vec3.length()) << std::endl;
//        }
//    }

    //float average_number_of_contacts = 4.1;
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

