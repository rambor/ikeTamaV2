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
/**
 * Random add/remove to generate initial model for seeded modeling
 * Only one PDB seed Allowed
 * Try to find the minimial set of lattice points that agree with atomistic P(r)
 * Uses Number of Contacts per bead in objective function, set eta to 0 for none
 *
 */
bool Anneal::createSeedFromPDB(PointSetModel *pModel, PofRData *pData, std::string name, std::string PDBFilename){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    contactCutOff = interconnectivityCutOff;
    auto lowTempStop = (double)highTempStartForCooling;

    // convert distances within the large search space to bins based on input P(R)-DATA file
    this->fillPrBinsAndAssignTotalBin( pModel, pData); // maxbin is set by diameter of the Universe
    unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    auto trueModelBeginIt = pModel->getSeedBegin();
    unsigned int workingLimit = pModel->getTotalInSeed();
    const unsigned int totalInModel = workingLimit;

    // copy trueModel into BeadIndices
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    /*
     * Create set of Anchors if exists
     * Anchors must be included in the model and connected to the model
     * Can have more than one anchor
     *
     * totalComponents default is 0, set by setAnchorPoints
     */
    std::set<unsigned int> excludeAnchorsList;
    for(unsigned int i=0; i < totalComponents; i++) {
        // the selected set of beads that make up each component will be held by Component object
        // add randomly selected indices to Component
        Component * component = &(components[i]);
        if (!component->anchorsEmpty()){
            for(auto it = component->getAnchors()->begin(); it != component->getAnchors()->end(); ++it) {  // find the anchor
                excludeAnchorsList.insert(*it);
            }
            component->writeAnchorsToFile("anchors_"+std::to_string(i+1));
        }
    }

    std::cout << " ESTIMATING MODEL FROM INPUT PDB => " << PDBFilename << std::endl;
    logger("TOTAL EXP N_S BINS", std::to_string(totalBins));
    logger("MAX MODEL N_S BINS", std::to_string(maxbin));
    logger("BINWIDTH", formatNumber(pData->getBinWidth(),2));
    logger("LATTICE SPACING", formatNumber(pModel->getBeadRadius()*2, 2));

    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    auto minWorkingLimit = (unsigned int)(0.17*workingLimit);

    std::vector<unsigned int> bead_indices(workingLimit); // large vector ~1000's
    std::vector<unsigned int> lowest_bead_indices(workingLimit); // large vector ~1000's
    std::vector<unsigned int> backUpState(workingLimit);

    // c-style is slightly faster for large vector sizes
    const unsigned int num = workingLimit;
    unsigned int * ptr = (num != 0) ? &bead_indices.front() : nullptr;
    for(unsigned int i = 0; i < num; i++) {
        ptr[i] = *(trueModelBeginIt+i);
    }

    // prepare bead_indices by copying in truemodel
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    // check for anchors, if anchors are present, move them out of the core model (anchors are only considered in the component)
//    std::cout << "                       ANCHORS ?  => " << excludeAnchorsList.size() << endl;
//    std::cout << "              INITIAL WORKINGLIMIT   " << workingLimit << endl;
//    if(excludeAnchorsList.size() > 0){
//        for(int i=0; i<workingLimit; i++){
//            if (excludeAnchorsList.find(bead_indices[i]) != excludeAnchorsList.end()){
//                // move to workingLimit and decrement
//                int decrement=1;
//                while( !( excludeAnchorsList.find(*(bead_indices.begin()+(workingLimit-decrement))) == excludeAnchorsList.end() ) ){
//                    decrement--;
//                }
//                cout << " FOUND ANCHOR AND SWAPPING TO END " << bead_indices[i] << endl;
//                std::iter_swap(bead_indices.begin() + i, bead_indices.begin() + (workingLimit-decrement));
//                workingLimit -= decrement;
//            }
//        }
//    }
//    std::cout << " (AFTER ANCHOR CHECK) WORKINGLIMIT  " << workingLimit << endl;
    const unsigned int deadLimit = workingLimit;; // as bead indices are discarded, set upper limit of vector
    //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // Surface area calculations
    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
    std::vector<float> weights(workingLimit);

    //float deltar = pModel->getBeadRadius()*0.74f;
    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    std::vector<Eigen::Vector3f> coordinates(workingLimit);

    for(unsigned int i=0; i < workingLimit; i++){
        const vector3 & pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
        unsorted_bead_indices[i] = bead_indices[i];
    }

    const auto sa = (float)(pModel->getBeadRadius()*pModel->getBeadRadius()*4.0*M_PI);
    float startingStoV = surfaceToVolume(workingLimit, weights, coordinates);

    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.end());
    double inv_kb_temp = 1.0/(double)highTempStartForCooling;

    auto lowerN = (unsigned int)(0.1*workingLimit);
    unsigned int upperN = workingLimit;

    logger("TOTAL LATTICE IN SEED", std::to_string(workingLimit));
    logger("LATTICE LIMITS", std::to_string(lowerN));

    // randomize and take the workingLength as first set, shuffling takes about 10x longer than copy and sort

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isConnected = true;
    // bead_indices contains only the indices that relate to the input PDB
    logger("CHECKING CONNECTIVITY", "");
    logger("NUMBER OF COMPONENTS", std::to_string(currentNumberOfComponents));

    //const int components = tempNumberOfComponents;
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    auto beginBinCount = binCount.begin();
    auto endBinCount = binCount.end();

    // calculate KL against known
    // calculateKLEnergy populates binCount based on model in bead_indices
    // create custom D_KL function supplying Pr_of_PDB as target
    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel);
    float invShannonNumber = 1.0f/(float)pData->getShannonBins();
    float testKL, currentKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB)*invShannonNumber;

    //float tempTotalContactEnergy = calculateTotalContactEnergy(&bead_indices, workingLimit, pModel, pDistance);
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    unsigned int lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

    // int priorWorkingLimit;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    //output for plotting
    bool isUpdated = false;

    auto coupons = (unsigned int)(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5d);
    unsigned int updateCount = ccmultiple*coupons;
    float seedHighTempRounds = updateCount/0.65f;///0.1;
    unsigned int updatedCount = 0;

    //unsigned int seedHighTempRounds = 31*deadLimit;
    float startContactSum=0.0f;
    for (unsigned int i=0; i<lowestWorkingLimit; i++){
        startContactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }

//    double lowestRunningContactSum, runningContactsSum = calculateTotalContactSumPotential(&beads_in_use_tree, pModel);
//    double etaConstant = eta*currentKL/((1.0-eta)*runningContactsSum);
//    double delTotalContactEnergy, totalContactEnergy = etaConstant*runningContactsSum;
    //double etaFactor = 1.0/std::pow(10000, 1.0/(seedHighTempRounds/(float)deadUpdate));
    startContactSum = startContactSum/(float)workingLimit;

    float this_energy, startKL = currentKL, lowestKL = currentKL;
    float current_energy = currentKL + (float)(lambda*connectivityPotential(currentNumberOfComponents));
    char addRemoveText[50];

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int counter=1, original;
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    unsigned int high = 0;

    for (; high < seedHighTempRounds; high++){ // iterations during the high temp search

        std::cout << "******************************************************************************" << std::endl;

        if (distribution(gen) >= 0.51 && workingLimit < (deadLimit-3)){ // ADD BEAD?

            std::cout << "*******************                  ADD                   *******************" << std::endl;
            std::sprintf(addRemoveText, " ");

            original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            auto itIndex = std::find(bead_indices.begin(), bead_indices.end(), original);
            auto distance = (unsigned int) std::distance(bead_indices.begin(), itIndex);

            while(original == noNeigborIndex || distance >= deadLimit){
                original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                itIndex = std::find(bead_indices.begin(), bead_indices.end(), original);
                distance = (unsigned int)std::distance(bead_indices.begin(), itIndex);
            }

            // select first element from randomized active_set
            // check if available neighbor can be added
            // make the swap at the border (workingLimit)
            addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);
            addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

            testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB)*invShannonNumber;

            // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
            //this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);
            this_energy = testKL;
            beads_in_use_tree.insert(original);

            if ( this_energy < current_energy || ( exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                currentNumberOfComponents = 1;
                isUpdated = true;
                eulerTour.addNode(original, pModel);
                std::sprintf(addRemoveText, "     ADD => %i", 1);
            } else { // undo changes (rejecting)
                beads_in_use_tree.erase(original);
                restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                // restore surface area changes
//                    unsorted_bead_indices.pop_back();
//                    weights.pop_back();
//                    coordinates.pop_back();
            }

        } else { // REMOVE BEADS?
            std::cout << "*******************                 REMOVE                 *******************" << std::endl;
            // test for deletion
            std::sprintf(addRemoveText, "     REMOVE => %i", 1);
            original = bead_indices[randomIndex(gen)];

            bool tourtest = true;
            while(tourtest){
                if (eulerTour.removeNode(original) == 1){
                    tourtest = false;
                } else {
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[ randomIndex(gen) ]; // potential to select the same thing twice
                }
            }

            if (!tourtest){
                // grab from randomized active_indices list
                removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);

                testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB)*invShannonNumber;
                this_energy = testKL;

                beads_in_use_tree.erase(original);

                if ( this_energy < current_energy || ( exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    currentNumberOfComponents = 1;
                    isUpdated = true;
                } else { // undo changes and move to next bead (rejecting)
                    auto beginIt = bead_indices.begin();
                    beads_in_use_tree.insert(original);
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.addNode(original, pModel);
                }
            } else {
                eulerTour.addNode(original, pModel);
            }
        }



        if (isUpdated){
            currentKL = testKL;
            current_energy = this_energy;

            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

        } else {
            acceptRate = inv500*(499*acceptRate);
        }

        updateASATemp(high, seedHighTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        std::cout << "*******************                                        *******************" << std::endl;
        printf("      TEMP : %-.3E => INVKBT %-.3E ACCEPTRATE %.2f\n", lowTempStop, inv_kb_temp, acceptRate);
        printf("  MAXSTEPS : %.0f ( %i ) \n", seedHighTempRounds, high);
        printf("     GRAPH : %i UPDATED %i\n", currentNumberOfComponents, updatedCount);
        printf("     LIMIT : %5i ( >= MIN: %i %i)  \n", workingLimit, minWorkingLimit, totalInModel);
        printf("      D_KL : %.4E SCORE : %.4E \n", currentKL, current_energy);


        if (currentKL < lowestKL){
            lowestKL = currentKL;
            lowestWorkingLimit = workingLimit;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            updatedCount++;
        }

        counter++;
    } // end of HIGH TEMP EQUILIBRATION

    /*
     * recalculate Surface and Volume of final model
     */
    weights.resize(workingLimit);
    coordinates.resize(workingLimit);

    for(unsigned int i=0; i < workingLimit; i++){
        const vector3 & pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
    }
    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    float finalStoV = surfaceToVolume(workingLimit, weights, coordinates);


    // average number of contacts per bead
    float contactSum=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        contactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }


    // determine the distribution of contacts
    // distribution matching
    unsigned int totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    logger("DISTRIBUTION OF CONTACTS", "-- CONTACTS --");
    for (unsigned int c=1; c<13; c++){
        unsigned int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        logger(std::to_string(c), formatNumber((totalContactsAt/(double)totalCounts), 2));
    }


    float average_number_of_contacts = contactSum/(float)workingLimit;

    logger("KL DIVERGENCE", "");
    logger("INITIAL D_KL", formatNumber(startKL,7));
    logger("LOWEST D_KL", formatNumber(lowestKL,7));
    logger("", "AVERAGE CONTACTS");
    logger("INITIAL", formatNumber(startContactSum,2));
    logger("FINAL (bin[one])", formatNumber(((2.0*binCount[0])/(double)workingLimit),2));
    // this is fixed model for initial high temp search?
    // set this as the seed
    //    pModel->setBeadAverageAndStdev(average_x, stdev);
    //    pModel->setReducedSeed(workingLimit, bead_indices);
    //    pModel->writeModelToFile(workingLimit, bead_indices, name);
    //    pModel->setBeadAverageAndStdev(average_x, stdev);

    // remap to bead universe for storage and printing
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    bead_indices.resize(totalBeadsInSphere);
    ptr = &bead_indices.front();
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::vector<unsigned int>::iterator it;
    for (unsigned int i=0; i<workingLimit; i++){
        it = std::find(bead_indices.begin(), bead_indices.end(), lowest_bead_indices[i]); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(bead_indices.begin()+i, it);
    }

    //pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);

    logger("LOWEST WORKINGLIMIT", std::to_string(lowestWorkingLimit));
    // set seed model to be invariant during reconstruction

    char flags[25];
    std::sprintf(flags, "qhull s FA");
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    /*
     * REMOVE ALL ANCHORS FROM REDUCED MODEL
     * Anchors are not part of the core reduced model but will be part of the model that is built ab initio
     */
    std::cout << "                       ANCHORS ?  => " << excludeAnchorsList.size() << std::endl;
    if(!excludeAnchorsList.empty()){
        for (auto eit : excludeAnchorsList){
            auto found = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, eit);
            auto dis = (unsigned int)std::distance(bead_indices.begin(), found);
            if (dis < workingLimit){ // move it
                std::cout << *it <<  " FOUND ANCHOR AND SWAPPING TO END " << " " << workingLimit << std::endl;
                std::iter_swap(found, bead_indices.begin() + workingLimit-1);
                workingLimit--;
            }
        }
    }

    logger("FINAL WORKINGLIMIT", std::to_string(workingLimit));
    logger("INITIAL S-to-V", formatNumber(startingStoV,2));
    logger("FINAL S-to-V", formatNumber(finalStoV,2));

    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    pModel->setReducedSeed(workingLimit, bead_indices);

    std::string nameOfModel = pModel->writeModelToFileBare(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            name,
            this,
            high,
            cvx,
            average_number_of_contacts);


    //create a distribution from the reduced model
    std::vector<unsigned int> contactsDistributionSeed(13);
    std::fill(contactsDistributionSeed.begin(), contactsDistributionSeed.end(), 0);

    for(unsigned int i = 0;i<workingLimit; i++){
        contactsDistributionSeed[numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i])]++;
    }

    // update distribution
    setContactsDistribution(contactsDistributionSeed);

    return true;
}


/**
 * Refinement will take the best fitting model from initialModelCVXHullSeeded and refine using Adaptive Simulated Annealing
 * The routine only modifies lattice positions defined by each instance of Component.
 *
 */
std::string Anneal::refineHomogenousBodyASACVXSeeded(PointSetModel *pModel, PofRData *pData, std::string outputname){
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;

    auto lowTempStop = (double)highTempStartForCooling;
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    // copy Starting_Set from initial model
    pModel->copyStartingModelIntoVector(bead_indices);
    std::set<unsigned int> reduced_model_tree(pModel->getReducedSeedBegin(), pModel->getReducedSeedEnd()); // reducedModel
    std::vector<unsigned int> reduced_model_vec(reduced_model_tree.size());
    std::copy(reduced_model_tree.begin(), reduced_model_tree.end(), reduced_model_vec.begin());

    std::set<unsigned int> seed_model;
    pModel->copySeedVectorIntoSet(seed_model);

    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    // reset components
    std::set<unsigned int> anchorsInUse;
    for (auto &it : components){
        it.copyBestToInUse();
        if (!it.anchorsEmpty()){
            for (auto sit : *it.getCenteredAnchors()){
                anchorsInUse.insert(sit);
            }
        }
    }

    for(unsigned int c=0; c<totalComponents; c++){
        Component * comp = &(components[c]);
        auto pSet = components[c].getBeadsInUse();
        for(auto const &i : *pSet){ // iterate over the set

            auto locale = (unsigned int)std::distance(bead_indices.begin(), std::find(bead_indices.begin(), bead_indices.end(), i) );
            if (locale >= workingLimit){
                std::cout << comp->getID() << " | " << ". pSET Bead not found within WorkingLimit " << i << " "  << locale << " ( " << workingLimit << " )" << std::endl;
                std::cout << "       isAnchor " << comp->isAnchor(i) << std::endl;
                std::cout << " centeredAnchor " << comp->isCenteredAnchor(i) << std::endl;
                std::cout << " workingLimit " << reduced_model_tree.size() << " + " << pSet->size() << " = " << workingLimit << std::endl;
                throw std::invalid_argument( "COMPONENT INDICES NOT FOUND : CRITICAL ERROR\n");
            }
        }
        //comp->printCenteredAnchors();
    }

    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    this->fillPrBinsAndAssignTotalBin(pModel, pData);
    unsigned short int * const pBin = pModel->getPointerToBins(); // pBin is required for pre-computed spaces
    /*
     * create cross-validation set or working observed probability distribution that encompasses search space
     */
    pData->setScoringFunction(maxbin);
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::cout << "    TOTAL EXP N_S BINS: " << totalBins << std::endl;
    std::cout << "    MAX MODEL N_S BINS: " << maxbin << std::endl;
    std::cout << "              BINWIDTH: " << pData->getBinWidth() << std::endl; // lattice should be limited by binwidth

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel);
    float testKL = 0.0f, currentKL = (float)pData->getScore(binCount);

    // Surface area calculations
    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
    std::vector<float> weights(workingLimit);

    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    std::vector<Eigen::Vector3f> coordinates(workingLimit);

    for(unsigned int i=0; i<workingLimit; i++){
        vector3 pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
        unsorted_bead_indices[i] = bead_indices[i];
    }

    const auto sa = (float)(4.0f*M_PI*(pModel->getBeadRadius() + delta_r)*(pModel->getBeadRadius() + delta_r));
    float tempStoV, currentStoV = surfaceToVolume(workingLimit, weights, coordinates);

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    std::uniform_real_distribution<float> distribution(0.0,1.0);

    float tempAverageContacts;

    // coupon collector's problem
    auto coupons = (unsigned int)(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5d);
    unsigned int updateCount = ccmultiple*coupons;

    float step_limit = (updateCount < 10000) ? 27000 : updateCount;

    std::vector<float> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    float * pTempDuringRun = &tempDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();
    unsigned int tempNumberOfComponents = currentNumberOfComponents;

    bool isUpdated = false;
    float acceptRate = 0.5, inv500 = 1.0f/500.0f;

    double inv_kb_temp = 1.0/lowTempStop;
    float this_energy, lowestKL = currentKL;
    char addRemoveText[50];

    // calculate initial/current energy of the system
    float componentPotential=0;
    for(int i=0; i < totalComponents; i++) {
        componentPotential+=components[i].potential();
    }

    // want D_KL to settle < 10^-5 to 10^-6
    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);
    const auto eightyfivePercent = (unsigned int)(0.85*step_limit);


    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;

    float muConstant = mu*0.1f;//0.0001;
    float volume_energy, test_volume, temp_component_volume_energy, current_component_volume_energy, currentComponentPotential = componentPotential;
    float tempConnectivityPotential, currentConnectivityPotential = connectivityPotentialPhases(currentNumberOfComponents);
    float current_energy = currentKL + lambda*currentConnectivityPotential + beta*currentComponentPotential;// + sToVConstant*currentStoV;
    float lowest_energy = current_energy;


    /*
     * check Component set that each node is present in bead_indices
     */
    for(auto &comp : components){
        auto pSet = comp.getBeadsInUse();
        auto pHull = comp.getHull();

        for(auto i : *pSet){
            if (beads_in_use_tree.find(i) == beads_in_use_tree.end()){
                std::cout << " Beads not found " << i << std::endl;
                std::cout << " pSet " << pSet->size() << " " << beads_in_use_tree.size() << std::endl;
                throw std::runtime_error("refineHomogenousBodyASACVXSeeded: COMPONENT INDICES NOT FOUND IN SET : CRITICAL ERROR\n");
            }

            auto locale = (unsigned int)std::distance(bead_indices.begin(), std::find(bead_indices.begin(), bead_indices.end(), i));
            if (locale >= workingLimit){
                std::cout << " Beads not found in vector " << i << std::endl;
                std::cout << " pSet " << pSet->size() << " " << beads_in_use_tree.size() << std::endl;
                throw std::runtime_error("refineHomogenousBodyASACVXSeeded: COMPONENT INDICES NOT FOUND WITHIN WL : CRITICAL ERROR\n");
            }
        }

        /*
         * for each component, set : hull, contacts distribution and divergence
         */
        populateLayeredDeadlimitUsingSet(beads_in_use_tree, *pHull, pModel);
        populateContactsDistribution(*comp.getContactsDistribution(), pSet, pModel);
        comp.setContactsKLDivergence(calculateKLDivergenceContactsDistribution(*comp.getContactsDistribution()));

        if (comp.getCurrentCVXVolume() <= 0){
            comp.setCurrentCVXVolume();
        }
    }

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

//    std::normal_distribution<float> volumeGen(pModel->getVolumeAverage(), pModel->getVolumeStdev());

    std::clock_t startTime;
    int attempts=0, failures=0;

    double runtime, tempKLDivContacts;
    unsigned int numberOfCoolingTempSteps, componentIndex, neighborLimit = pModel->getNeighborLimit();

    std::uniform_int_distribution<unsigned int> compIndex(0, totalComponents-1); // guaranteed unbiased

    for(numberOfCoolingTempSteps = 0; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        std::sprintf(addRemoveText, "");
        /*
         * randomly select a component to modify
         * calculate current volume and contacts energy
         */
        componentIndex = (totalComponents == 1) ? 0 : compIndex(gen);
        Component * pComponent = &components[componentIndex];
        std::set<unsigned int> * pComponentSetTree = pComponent->getBeadsInUse();

        current_component_volume_energy = muConstant*pComponent->calculateVolumeEnergy(pComponent->getCurrentCVXVolume());

        double contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*pComponent->getContactsKLDivergence());
        double currentKLDivContactsEnergy = contactDistributionWeight*pComponent->getContactsKLDivergence();

        std::vector<unsigned int> selections(pComponentSetTree->begin(), pComponentSetTree->end());
        std::shuffle(selections.begin(), selections.end(), gen);
        auto setLimit = (unsigned int)pComponentSetTree->size();

        // set the hull
        auto pHull = pComponent->getHull();
        auto hullsize = (unsigned int)pHull->size();
        std::uniform_int_distribution<unsigned int> randomHull(0,hullsize-1); // guaranteed unbiased

        std::vector<double> tempContactsDistributionOfModel(pComponent->getContactsDistribution()->begin(), pComponent->getContactsDistribution()->end());

        startTime = std::clock();

        if ( distribution(gen) < percentAddRemove ) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;

            if ( distribution(gen) < 0.63 ) {
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                /*
                 * pick random bead to add to treats each bead equally which means that the neighbor we pick is biased
                 * to denser regions.
                 * picking random bead within hull that is not in use biases to non-selected regions
                 */
                unsigned int testIndex = selections[0];
                for(auto & sel : selections){
                    testIndex = getUseableNeighborWeighted(&beads_in_use_tree, pModel, sel);
                    if (testIndex != pModel->getNeighborLimit()){
                        break;
                    }
                }

//                auto set_it = pHull->begin();
//                std::advance(set_it, randomHull(gen));
//                unsigned int testIndex = *set_it;
                auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), testIndex);

                //auto itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);

                addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);
                addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                testKL = (float)pData->getScore(binCount);

                pComponent->addLatticePoint(testIndex);
                addToContactsDistribution(testIndex, tempContactsDistributionOfModel, pComponentSetTree, pModel);

                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
                double tempKLDivContactsEnergy = contactDistributionWeight*tempKLDivContacts;

                test_volume = pComponent->calculateCVXVolume();
                temp_component_volume_energy = muConstant*pComponent->calculateVolumeEnergy(test_volume);

                // may be unnecessary
                componentPotential=0;
                for(unsigned int i=0; i < totalComponents; i++)
                    componentPotential+=components[i].potential();

                // determine Euler tours
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                tempConnectivityPotential = connectivityPotentialPhases(tempNumberOfComponents);

                this_energy = testKL +
                        lambda*tempConnectivityPotential +
                        beta*componentPotential +
                        temp_component_volume_energy +
                        tempKLDivContactsEnergy;

                if (this_energy  < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                    beads_in_use_tree.insert(testIndex);
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     ADD => %i", testIndex);
                } else { // undo changes (rejecting)
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                    pComponent->removeLatticePoint(testIndex);
                    std::sprintf(addRemoveText, "     ADD => FAILED");
                }

            } else { // REMOVE BEADS?
                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                startTime = std::clock();
                // go through entire list
                unsigned int testIndex;
                unsigned int currentNumberOfToursInComponent = pComponent->getTotalNumberOfComponents();
                bool removeit = false;

                for(auto & select : selections){
                    if( (reduced_model_tree.find(select) == reduced_model_tree.end()) && canRemoveIfNotCenteredAnchor(select)){

                        if (pComponent->removeLatticePoint(select) <= currentNumberOfToursInComponent){
                            tempNumberOfComponents = eulerTour.removeNode(select);
                            if (tempNumberOfComponents <= currentNumberOfComponents){
                                removeit = true;
                                testIndex = select;
                                break;
                            }
                            eulerTour.addNode(select, pModel);
                        }
                        pComponent->addLatticePoint(select);
                    }
                }


                if (removeit){ // may not find a suitable node to remove, so must abort if no node found
                    removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);

                    testKL = (float)pData->getScore(binCount);

                    componentPotential=0;
                    for(unsigned int i=0; i < totalComponents; i++) { // the selected set of beads that make up each component will be held by Component object
                        componentPotential+=components[i].potential();
                    }


                    removeFromContactsDistribution(testIndex, tempContactsDistributionOfModel, pComponentSetTree, pModel);
                    tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
                    double tempKLDivContactsEnergy = contactDistributionWeight*tempKLDivContacts;

                    tempConnectivityPotential = connectivityPotentialPhases(tempNumberOfComponents);

                    test_volume = pComponent->calculateCVXVolume();
                    temp_component_volume_energy = muConstant*pComponent->calculateVolumeEnergy(test_volume);

                    this_energy = testKL +
                                  lambda*tempConnectivityPotential +
                                  beta*componentPotential +
                                  temp_component_volume_energy +
                                  tempKLDivContactsEnergy;

                    if (this_energy < current_energy || exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)) {
                        beads_in_use_tree.erase(testIndex);
                        std::sprintf(addRemoveText, "     REMOVED => %i", testIndex);
                        isUpdated = true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        pComponent->addLatticePoint(testIndex);

                        auto beginIt = bead_indices.begin();
                        auto beginBinCount = binCount.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                        std::sprintf(addRemoveText, "     REMOVED => FAILED");
                    }
                }
            }
        } else {
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            startTime = std::clock();

            unsigned int swap1 = selections[0];
            unsigned int currentNumberOfToursInComponent = pComponent->getTotalNumberOfComponents();

            bool lost = false;

//            for(const auto & select : selections){
//                swap1 = select;
//                if (!pComponent->isCenteredAnchor(swap1)){
//                    if ( pComponent->removeLatticePoint(swap1) <= currentNumberOfToursInComponent ){
//                        if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
//                            lost = true;
//                            break;
//                        }
//                        eulerTour.addNode(swap1, pModel);
//                    }
//                    pComponent->addLatticePoint(swap1);
//                }
//            }

            for(const auto & select : selections){
                swap1 = select;
                if (((1-getProbabilityOfSelectedIndex(&beads_in_use_tree, pModel, swap1)) < distribution(gen)) && !pComponent->isCenteredAnchor(swap1)){
                    if ( pComponent->removeLatticePoint(swap1) <= currentNumberOfToursInComponent ){
                        if (eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            lost = true;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                    pComponent->addLatticePoint(swap1);
                }
            }


            // the component's beads_in_use tree has been modified
            if (lost){ // if I can't find a suitable bead to move, skip
                auto itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
                // remove selected index from P(r)
                removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                // find better position
                beads_in_use_tree.erase(swap1);
                /*
                 * random grab neighbor from selections
                 */
                unsigned int neighbor = neighborLimit;
                for(unsigned int i=(setLimit-1); i > 0; --i){ // work backwards?
                    unsigned int newPosition = selections[i];
//                    if (newPosition != swap1 ){
//                        neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, newPosition);
//                        if (neighbor != neighborLimit && neighbor != swap1 && (seed_model.find(neighbor) == seed_model.end()) ){
//                            break;
//                        }
//
//                        if (neighbor != neighborLimit && neighbor != swap1  ){
//                            break;
//                        }
//                    }
                    neighbor = getUseableNeighborWeighted(&beads_in_use_tree, pModel, newPosition);
                    if (neighbor != neighborLimit && neighbor != swap1 && (seed_model.find(neighbor) == seed_model.end())){
                        break;
                    }
                }


                // make the swap, sort and update P(r)
                auto pSwap2 = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), neighbor);
                std::iter_swap(itIndex, pSwap2);
                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

                addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
                // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                testKL = (float)pData->getScore(binCount);
                // Euler Tours and Connectivity
                tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);
                pComponent->addLatticePoint(neighbor);

                populateContactsDistribution(tempContactsDistributionOfModel, pComponentSetTree, pModel);
                double tempKLDivContactsEnergy = contactDistributionWeight*calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

                componentPotential=0;
                for(unsigned int i=0; i < totalComponents; i++)  // the selected set of beads that make up each component will be held by Component object
                    componentPotential+=components[i].potential();

                // update volume
                test_volume = pComponent->calculateCVXVolume();
                temp_component_volume_energy = muConstant*pComponent->calculateVolumeEnergy(test_volume);

                tempConnectivityPotential = connectivityPotentialPhases(tempNumberOfComponents);

                this_energy = testKL +
                              lambda*tempConnectivityPotential +
                              beta*componentPotential +
                              temp_component_volume_energy +
                              tempKLDivContactsEnergy;

                if (this_energy < current_energy || exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)) {
                    beads_in_use_tree.insert(neighbor);
                    std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, neighbor);
                    isUpdated = true;
                } else {
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    // tour of the entire space
                    eulerTour.removeNode(neighbor);
                    eulerTour.addNode(swap1, pModel);
                    // tour of the component
                    pComponent->removeLatticePoint(neighbor);
                    pComponent->addLatticePoint(swap1);
                    beads_in_use_tree.insert(swap1);
                    std::sprintf(addRemoveText, "     SWAPPED => FAILED");
                }
            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        /*
         * TESTS
         */
        // failures usually mean backups were not properly placed, or indices were mis-added/removed in incorrect order
//        calculateModelPrDistribution(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        float testKL1 = pData->getScore(testBinCount);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)  ){
//            std::cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << std::endl;
//            return "stopped";
//        }
//        if (connectivityPotentialPhases(eulerTour.getNumberOfComponents()) > 0){
//            cout << " STOPPED TOO MANY COMPONENTS " << currentNumberOfComponents << " " << eulerTour.getNumberOfComponents() << endl;
//            string name = "failed_" + std::to_string(numberOfCoolingTempSteps);
//            //pModel->writeModelToFile(workingLimit, bead_indices, name);
//            return "Stopped";
//        }

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            currentConnectivityPotential=tempConnectivityPotential;
            currentComponentPotential = componentPotential;

            pComponent->setCurrentCVXVolume(test_volume);
            pComponent->setContactsKLDivergence(tempKLDivContacts);
            populateLayeredDeadlimitUsingSet(*pComponentSetTree, *pHull, pModel);

            current_component_volume_energy = temp_component_volume_energy;

            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;
            failures=0;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

            if (currentKL < lowestKL)
                lowestKL = currentKL;

            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), pComponent->getContactsDistribution()->begin());

        } else {
            acceptRate = inv500*(499*acceptRate);
            failures++;
        }

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        printf("      TEMP : %-.4E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %.5f FAILURES => %i  TIME : %.4f\n", acceptRate, failures, runtime);
        printf("     LIMIT : %5i  GRAPH => %3i (POT %.5F) \n", workingLimit, currentNumberOfComponents, currentConnectivityPotential);
        printf("    %s\n", addRemoveText);
        printf(" COMPONENT : %d VOL_E => %.3E VOL => %.3E %.3E \n", componentIndex, current_component_volume_energy, pComponent->getCurrentCVXVolume(), currentKLDivContactsEnergy);
        printf("     D_KL => %-4.3E ( %.4E ) ENRGY : %.4E \n", currentKL, lowestKL, current_energy);

        pTempDuringRun[numberOfCoolingTempSteps] = lowTempStop;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;


        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        //update for running average
        if (numberOfCoolingTempSteps < fifteenPercent){
            currentCDW *= alphaConstant;
        }
        // rescale etaFactor in small steps until final value
    } // end of steps

    // At end of each temp, update a probability model for volume?  Use this to select
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    printf(" NUMBER OF STEPS %i\n", numberOfCoolingTempSteps);
    printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    std::cout << "------------------------------------------------------------------------------" << std::endl;

    tempAverageContacts=0.0;
    for (int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << std::endl;
    std::cout << " POSITIONAL SUCCESS/FAILURE " << ((double)failures/(double)attempts)*100 << std::endl;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    // perform positional refinement until delta E stabilizes?
    int totalLessThanAverage=0, numberContacts;
    for (int i=0; i<workingLimit; i++){
        numberContacts = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        if (numberContacts < average_number_of_contacts ){
            totalLessThanAverage +=1;
        }
    }

    std::cout << " TOTAL LATTICE POINTS LESS THAN AVERAGE " << totalLessThanAverage << std::endl;

    for(int i=0; i < totalComponents; i++) {
        components[i].printConstraints(); // constraints the number of beads
    }
    //this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);

    //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);

    // if multithread, must put a lock on these two step
    char flags[25];
    sprintf(flags, "qhull s FA");
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    //pModel->setCVXHullVolume(calculateCVXHULLVolume(flags, &bead_indices, workingLimit, points, pModel));

    pData->printKLDivergence(binCount);
    // write components

    for (auto & pComp : components){
        pComp.writeToFile("component_" + pComp.getID() + "_"+outputname);
        std::cout << pComp.getID() << " TARGET NUMBER OF POINTS => " << pComp.getTargetNumberOfLatticePoints() << " ( " << pComp.getBeadsInUse()->size() << " )" << std::endl;
        pComp.printContactsPerCenteredAnchor();
        pComp.writeCenteredAnchorsToFile();
    }

    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            outputname,
            this,
            pData,
            numberOfCoolingTempSteps,
            cvx,
            average_number_of_contacts);


    return nameOfModel;
}

/**
 * Test if index is removable, returns true if index is not an anchor point
 * @param index
 * @return
 */
bool Anneal::canRemoveIfNotCenteredAnchor(unsigned int index) {

    for(unsigned int i=0; i < totalComponents; i++) {
        // the selected set of beads that make up each component will be held by Component object
        Component * comp = &components[i];
        if (comp->inUse(index)){
            return !comp->isCenteredAnchor(index);
        }
    }
    // if it doesn't belong to a component, assume it is part of the seed model
    return true;
}

/**
 * calculates connectivity potential as a simple quadratic
 * @param mainConnectivity
 * @return
 */
float Anneal::connectivityPotentialPhases(unsigned int mainConnectivity){
    float currentConnectivityPotential = (mainConnectivity-1.0f)*(mainConnectivity-1.0f);
    float temp;

    for(auto & comp : components){
        temp = comp.getTotalNumberOfComponents() - 1.0f;
        currentConnectivityPotential += temp*temp;
    }

    return currentConnectivityPotential;
}