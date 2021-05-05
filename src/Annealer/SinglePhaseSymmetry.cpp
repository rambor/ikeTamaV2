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
 * create initial model of complete object
 * does not require pre-computed distances, this is a direct method meaning all pairwise distances are calcualted on the fly
 */
bool Anneal::createInitialModelSymmetry(PointSetModel *pModel, PofRData *pData) {

    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    auto lowTempStop =  (double)0.0001;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    contactCutOff = pModel->getNeighborCutOffLimit();
    violation_limit = pModel->getBeadRadius();
    //violation_limit = sqrt(3)*pModel->getBeadRadius();
    //violation_limit = (float)(2.0/3.0*(sqrt(3)*pModel->getBeadRadius()));

    maxbin = pModel->getMaxBin(pData);
    // adjust this part
    maxbin += 2;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> lowestbinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

    logger("TOTAL EXP N_S BINS", formatNumber(totalBins));
    logger("MAX MODEL N_S BINS", formatNumber(maxbin));
    logger("BINWIDTH", formatNumber((float)pData->getBinWidth(), 2));
    logger("BEAD RADIUS", formatNumber((float)pModel->getBeadRadius(), 2));

    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // as bead indices are discarded, set upper limit of vector
    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> test_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> lowest_subUnit_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    const unsigned int num = totalBeadsInSphere;
    unsigned int * const ptr = (num != 0) ? subUnit_indices.data() : nullptr;
    for(unsigned int i = 0; i < num; i++) {
        ptr[i] = i;
    }

//    printSymModel(&subUnit_indices, totalBeadsInSphere, pModel, pData);
//    std::set<unsigned int> universe(subUnit_indices.begin(), subUnit_indices.end());
//    std::string universename = "uni";
//    pModel->writeSetToFile(universe, universename);
//    universe.clear();

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.331);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    const double targetSubUnitVolume = (2*lowerV+upperV)/3/(float)pModel->getNumberOfSubUnits();
    const auto radius_larger = (float)(pModel->getBeadRadius()*std::sqrt(7.0)/2.0);
    const auto lowerN = (unsigned int)(std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/3;
    const auto upperN = (unsigned int)(std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()));

//    const auto lowerN = (unsigned int)(std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/3;
//    const auto upperN = (unsigned int)(std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/2;

    unsigned int subUnitWorkingLimit = (upperN + lowerN)/2;
    const unsigned int multipleWL = subUnitWorkingLimit*7;

    logger("Bead Search Limit => LOWER N", formatNumber(lowerN));
    logger("Bead Search Limit => UPPER N", formatNumber(upperN));
    logger("INITIAL MODEL WL", formatNumber(subUnitWorkingLimit));
    logger("SYMMETRY", pModel->getSymmetry());
    logger("TOTAL SUBUNITS", formatNumber(pModel->getNumberOfSubUnits()));

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    logger("", "CREATING INITIAL RANDOM MODEL");
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());

    std::set<unsigned int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::set<unsigned int> hull;  // sort hull points into bead_indices
    recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, hull, pModel);
    std::string switched = "HULL";
    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    coordT points[3*(2*upperN)];
    char flags[25];
    std::sprintf(flags, "qhull s FA");

    float subunit_test_volume, subunit_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
//    float test_volume, current_volume = calculateCVXVolumeSymmetry(&subUnit_indices, subUnitWorkingLimit, pModel);

    // calculate starting energy
    unsigned int totalViolations, temp_violations;
    std::clock_t start = std::clock();
    calculateModelPrDistributionSym(&subUnit_indices, &binCount, subUnitWorkingLimit, totalViolations, pModel, pData );

    // fill binCount for first time
    double currentKL = pData->getScore(binCount);

    double hlambda = pData->getIsPr() ? 0.01d : 10.0d;
//    double muConstant = pData->getIsPr() ? 0.00001d : 0.1d; // chi tends to start in 100's whereas DKL is 0.1;
    double muPercent = 0.731;
    double muConstant = currentKL*muPercent/((1.0 - muPercent)*subunit_volume/(double)subUnitWorkingLimit);
//    const float targetVolume = (0.5f*(lowerV+upperV));///(float)pModel->getNumberOfSubUnits();
//    const auto lowerSubUnitV = (unsigned int)((pData->getVolume()  - pData->getVolume() *0.21)/(float)pModel->getNumberOfSubUnits());
//    const auto upperSubUnitV = (unsigned int)((pData->getVolume()  + pData->getVolume() *0.05)/(float)pModel->getNumberOfSubUnits());

    //float diff = biCameralVolumePotential(lowerSubUnitV, upperSubUnitV, targetSubUnitVolume);

    //float invTarVolmuConstant = muConstant/(targetVolume*targetVolume);
    const double invTarVolSubUnitmuConstant = muConstant/targetSubUnitVolume;
    unsigned int failures=0;

    bool updateCVX = false;

    logger("THREAD TIME (sec)", formatNumber( (float)((std::clock() - start)/(double) CLOCKS_PER_SEC), 6) );
    // high temp search
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");

    double testKL, test_energy; // sets alpha as a constant during high temp eq

    unsigned int lowestWorkingLimit = subUnitWorkingLimit;
    unsigned int tempNumberOfComponents, currentNumberOfComponents;

    EulerTour eulerTour(subUnit_indices.begin(), subUnitWorkingLimit, pModel);
    currentNumberOfComponents = eulerTour.getNumberOfComponents();

    double current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) ;
    //double current_energy = currentKL + hlambda*connectivityPotential(currentNumberOfComponents);
    std::cout << "       INITIAL => ENERGY :  " << current_energy << std::endl;
    std::cout << "       INITIAL =>   D_KL : " <<  currentKL << std::endl;
    std::cout << "       INITIAL => VOLUME : " << subunit_volume << std::endl;
    std::cout << "              VIOLATIONS : " << totalViolations << std::endl;
    std::cout << "       INITIAL        WL : " << subUnitWorkingLimit << std::endl;

    char addRemoveText[50];
    float sum_x_squared=0;
    unsigned int volumeCount=0, testIndex;
    float volumeSum=0, workingLimitSum=0;

    // what is tolerable number of violations per subsunit?
    current_energy += beta * totalViolations;

    double lowest_energy = current_energy;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(binCount.begin(), binCount.end(), lowestbinCount.begin());  // unaltered P(r)
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup cop

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int high, hullsize = (unsigned int)hull.size();
    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased

    //double temp_subunit_volume_energy, current_subunit_volume_energy = std::abs(subunit_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
    double temp_subunit_volume_energy, current_subunit_volume_energy = muConstant*subunit_volume/(double)subUnitWorkingLimit;

    unsigned int addFailed = 0, removeFailed =0;
    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        if (distribution(gen) < 0.17 || currentNumberOfComponents > 5){
            // create points from workSet to determine HULL
            std::sprintf(addRemoveText, "POSITIONAL");
            unsigned int position = randomIndex(gen);
            unsigned int swap1 = subUnit_indices[position];
            auto itIndex = subUnit_indices.begin() + position;

            // setup parameters for hull
            if ( (distribution(gen) < 0.37129) || currentNumberOfComponents > 2 ){ // only move CVX points if true

                unsigned int * const pSelections = subUnit_indices.data(); // initialized as empty in Model class
                for (unsigned int i = 0; i < subUnitWorkingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                    active_indices[i] = pSelections[i];
                }

                // needs to be optimized
                qh_new_qhull(3, subUnitWorkingLimit, points, 0, flags, nullptr, nullptr);
                vertexT * vertices = qh vertex_list;
                auto totalV = (unsigned int)(qh num_vertices);

                // only move CVX hull points
                std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
                for (unsigned int v=0; v < totalV; v++) { //
                    indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                    vertices = vertices->next;
                }

                qh_freeqhull(true);

                std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
                swap1 = indices_to_check[randomVertices(gen)];
                itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit, swap1);
            } // recalculate

            // find bead to swap in active set
            removeFromPrSym(swap1, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData); // remove selected index from P(r)
            //beads_in_use_tree.erase(swap1);
            eulerTour.removeNode(swap1);
            // find position to swap to
            auto set_it=hull.begin();
            std::advance(set_it, randomHull(gen));
            unsigned int neighbor = *set_it;
            while (neighbor == swap1 || beads_in_use_tree.find(neighbor) != beads_in_use_tree.end()) { // find a new position
                set_it=hull.begin();
                std::advance(set_it, randomHull(gen));
                neighbor = *set_it;
            }

            // make the swap, sort and update P(r)
            auto pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.end(), neighbor);

            std::iter_swap(itIndex, pSwap2);

            temp_violations = addToPrSym(neighbor, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
            testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
            //temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
            temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;

            //test_energy = testKL + hlambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;
            test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) +  beta*temp_violations;//diff*invTarVolmuConstant;

            if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                beads_in_use_tree.erase(swap1);
                beads_in_use_tree.insert(neighbor);
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                //beads_in_use_tree.insert(swap1);
                eulerTour.removeNode(neighbor);
                eulerTour.addNode(swap1, pModel);
            }

        } else { // add or remove

            // hard limit on number of beads
            if ( (distribution(gen) < 0.51 && subUnitWorkingLimit > lowerN) ||  subUnitWorkingLimit > upperN  ) { // REMOVE beads from sorted list into useable range < deadLimit
                // randomly swap positions with end of workingLength, could remove CVX Hull Point
                std::sprintf(addRemoveText, "  REMOVE  "); // must only consider removing beads that do not increase number of components

                testIndex = subUnit_indices[randomIndex(gen)];
                while(true){
                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(testIndex, pModel);
                    testIndex = subUnit_indices[randomIndex(gen)];
                }


                removeLatticePositionToModelSym(subUnit_indices, binCount, &subUnitWorkingLimit, &testIndex, pModel, pData);
                temp_violations = getViolations(&subUnit_indices, subUnitWorkingLimit, pModel);

                testKL = pData->getScore(binCount);

                subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
                //temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
                temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;
//                tempNumberOfComponents = eulerTour.removeNode(testIndex);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations;
//                test_energy = testKL + hlambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;

                if ( ( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy - temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.erase(testIndex);
                    updateCVX=true;
                } else { // undo changes and move to next bead (rejecting)
                    eulerTour.addNode(testIndex, pModel);
                    std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                    subUnitWorkingLimit += 1;
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); //copy to bin count
                    removeFailed += 1;
                }

            } else { // ADD beads
                std::sprintf(addRemoveText, "    ADD   ");

//                std::set<unsigned int>::const_iterator set_it=hull.begin();
//                std::advance(set_it, randomHull(gen));
//                testIndex = *set_it;
//
//                while ( beads_in_use_tree.find(testIndex) != beads_in_use_tree.end()) { // find a new position
//                    set_it=hull.begin();
//                    std::advance(set_it, randomHull(gen));
//                    testIndex = *set_it;
//                }


                // if we select a random bead to add to, we will bias to, comparatively, denser regions
                unsigned int testIndex = getUseableNeighborWeighted(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                while(testIndex == pModel->getNeighborLimit()){
                    testIndex = getUseableNeighborWeighted(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
                }

                auto itIndex = std::find(subUnit_indices.begin()+ subUnitWorkingLimit, subUnit_indices.end(), testIndex);
                std::iter_swap(subUnit_indices.begin() + subUnitWorkingLimit, itIndex); // this swaps to working position and changes backup
                // increment workingLimit to include new position

                if (itIndex == subUnit_indices.end()){
                    std::cout << "Not found";
                    exit(0);
                }

                subUnitWorkingLimit += 1;

                temp_violations = addToPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);

                testKL = pData->getScore(binCount);

                subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
                //temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
                temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel); // want a single tour (so implies connectivity)

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations;// + diff*invTarVolmuConstant + beta*temp_violations;

                if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.insert(testIndex);
                    updateCVX=true;
                } else { // undo changes (rejecting)
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                    addFailed += 1;
                }
            }
        }



        if (updateCVX){

            current_energy = test_energy;
            currentKL = testKL;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            subunit_volume = subunit_test_volume;
            current_subunit_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            // randomly select points within convex hull or neighbor
            if ( hullsize < subUnitWorkingLimit/2 || currentNumberOfComponents == 1){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
                switched = "NEAREST NEIGHBOR";
            } else {
                switched = "HULL";
                recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, hull, pModel);
            }

            hullsize = hull.size();
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);

            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<unsigned int>(0,subUnitWorkingLimit-1); // guaranteed unbiased

            if (current_subunit_volume_energy > 1.5*currentKL && currentNumberOfComponents == 1){
                muConstant = currentKL*muPercent/((1.0 - muPercent)*subunit_volume/(double)subUnitWorkingLimit);
                current_subunit_volume_energy = muConstant*subunit_volume/(double)subUnitWorkingLimit;
            }
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASAConstantTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts

        std::printf("*******************             %s                 ******************* \n", addRemoveText);
        std::printf("      MAXSTEPS : %7i (%7i) ACCEPTRATE : %.3f TEMP : %.2E\n", highTempRounds, high, acceptRate, lowTempStop);
        std::printf("         GRAPH : %7i  MODE => %s \n", currentNumberOfComponents, switched.c_str());
        std::printf("         LIMIT : %7i UPPER <= %i  LOWER >= %i (%i)\n", subUnitWorkingLimit, upperN, lowerN, hullsize);
        std::printf("        VOLUME : %.0f (%.0f)  MU => %.1E  MU*VOL => %.2E\n", subunit_volume, targetSubUnitVolume, muConstant, current_subunit_volume_energy);
        std::printf("    VIOLATIONS : %7d  BETA => %.1E (%i) %i | %i\n", totalViolations, beta, volumeCount, addFailed, removeFailed);
        std::printf("          D_KL : %.4E ENRGY: %.4E (%.4E) \n", currentKL, current_energy, lowest_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy ){
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            std::copy(binCount.begin(), binCount.end(), lowestbinCount.begin());  // unaltered P(r)
            lowestWorkingLimit = subUnitWorkingLimit;

            lowest_energy = current_energy;
            workingLimitSum += subUnitWorkingLimit;
            sum_x_squared += subUnitWorkingLimit*subUnitWorkingLimit;
            volumeSum += subunit_volume;
            volumeCount++;
        }

        if (currentNumberOfComponents > 61 && high%multipleWL == 0){
            hlambda = currentKL*0.971/((1.0 - 0.971)*(currentNumberOfComponents-1)*(currentNumberOfComponents-1));
            current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + beta*totalViolations;
        }



//        checkSetAndVector(subUnitWorkingLimit, &subUnit_indices, &beads_in_use_tree);
//
//        if (checkForRepeats(subUnit_indices)){
//            std::cout << " STOPPED POSITIONAL " << " WL: " << subUnitWorkingLimit << " D_KL " << currentKL << " <=> " << std::endl;
//            exit(0);
//        }

//        violations = getViolations(&subUnit_indices, subUnitWorkingLimit, pModel, pData);
//        if (violations != totalViolations){
//            std::cout << " totalViolations " << totalViolations << " actual " << violations << std::endl;
//            testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, violations, pModel, pData );
//            if (currentKL != testKL){
//                printf("*******************             %s                 ******************* \n", addRemoveText);
//                std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//                std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//                printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//            }
//            exit(0);
//        }
//
//        calculateModelPrDistributionSym(&subUnit_indices, &testBinCount, subUnitWorkingLimit, temp_violations, pModel, pData );
//        testKL = pData->getScore(testBinCount); // 100x faster than calculateKLEnergySymmetry
//        if (currentKL != testKL){
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//            std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//            printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//            exit(0);
//        }
//
//        testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, temp_violations, pModel, pData );
//        if (currentKL != testKL){
//            // if this fails, update of P(r) is wrong or subUnit_indices is corrupted
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//            std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//            printf("CurrentKL %.8E  TestKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//
//            int bins = testBinCount.size();
//            std::cout << " BIN : " << std::endl;
//            for(int i=0; i<bins; i++){
//                if (testBinCount[i] != binCount[i]){
//                    std::cout << i << " : " << testBinCount[i] << " == " << binCount[i] << std::endl;
//                }
//            }
//
//            exit(0);
//        }

    } // end of HIGH TEMP EQUILIBRATION


    switched = "HULLFinal";
    recalculateDeadLimit(lowestWorkingLimit, lowest_subUnit_indices, hull, pModel);
    pModel->writeSetToFile(hull, switched);


    highTempStartForCooling = (float)lowTempStop;
    // pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage = workingLimitSum/(float)volumeCount;
    float volumeStdev = std::sqrt(sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage);

    // remove points close to hull
    std::string tempName = "initial_CVX_subunit_" + filenameprefix;

    if (currentNumberOfComponents > 1){
        pModel->setStartingSet(lowest_subUnit_indices);
        pModel->setStartingWorkingLimit(lowestWorkingLimit);
        pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, lowestWorkingLimit, lowest_subUnit_indices, lowestbinCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);
    } else {
        pModel->setStartingSet(subUnit_indices);
        pModel->setStartingWorkingLimit(subUnitWorkingLimit);
        pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, lowestbinCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);
    }

    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);

    temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel);
    logger("VIOLATIONS (FINAL)", std::to_string(temp_violations));
    pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, binCount, "final_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);


    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    logger("AVERAGE # BEADS HIGH TEMP SELECTION", std::to_string((int)volumeAverage));
    logger("SIGMA", std::to_string((int)volumeStdev));
    std::cout << "*******************                                        *******************" << std::endl;

    EulerTour finalEulerTour(lowest_subUnit_indices.begin(), lowestWorkingLimit, pModel);
    currentNumberOfComponents = finalEulerTour.getNumberOfComponents();

    if (currentNumberOfComponents == 1 ) {
        return true;
    } else {
        tempName = "failed_subunit_" + filenameprefix;
        pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, binCount, "failed_sym" + filenameprefix, this, pData, high, subunit_volume, 4.1);

        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return false;
    }
}



/**
 * Input model must be aligned to a subunit
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineSymModel(PointSetModel *pModel, PofRData *pData, std::string nameTo){
    const std::string sym = pModel->getSymmetry();
    bool once = true;
    const char * pnameTo = nameTo.c_str() ;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000000001 : (double)highTempStartForCooling;
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    /*
     * create set of beads to select from based on all those that are in contact with selected set
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    std::set<unsigned int> original_hull(beads_in_use_tree.begin(), beads_in_use_tree.end());

//    this->updateContactsDistribution(&beads_in_use_tree, pModel);

    if (maxbin < 2){
        throw std::invalid_argument(" INCORRECT MAXBIN NOT INITIALIZED PROPERLY");
    }

    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin);   // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);      // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", std::to_string(totalBins));
    logger("MAX MODEL N_S BINS", std::to_string(totalBins));
    logger("BINWIDTH (Angstroms)", formatNumber(pData->getBinWidth(), 3));

    unsigned int temp_violations, totalViolations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, totalViolations, pModel, pData );
    // fill binCount for first time
    double testKL, currentKL = pData->getScore(binCount);

    char flags[] = "qhull FA"; // CVX HULL STUFF

    // Surface area calculations
//    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
//    std::vector<float> weights(workingLimit);
//    probe_radius = pModel->getBeadRadius() + delta_r;
//    std::fill(weights.begin(), weights.end(), probe_radius);
//    std::vector<Eigen::Vector3f> coordinates(workingLimit);
//
//    for(unsigned int i=0; i<workingLimit; i++){
//        vector3 pBeadVec = pModel->getBead(bead_indices[i])->getVec();
//        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
//        unsorted_bead_indices[i] = bead_indices[i];
//    }

    // based on the subunit estimated from initial modeling, use as a upper limit of the volume
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float subunit_test_volume;
    const float target_volume = 0.8f*current_volume;
    double muConstant = 0.000001;
    const float invTarVolSubUnitmuConstant = 0;//muConstant/target_volume;
    double temp_subunit_volume_energy=0, current_volume_energy = 0;

    const float targetSubUnitVolume = (upperV)/(float)pModel->getNumberOfSubUnits();
    logger("TARGET VOLUME (A^3)", formatNumber(target_volume, 3));
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    float tempAverageContacts;
    // coupon collector's problem
    const unsigned int updateCount = ccmultiple*getCoupons(3*workingLimit);
    float step_limit = (updateCount < 51357) ? 51357 : (float)updateCount;

    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int tempNumberOfComponents = 0, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float acceptRate = 0.5;
    const float inv500slash499 = 499.0f/500.0f, inv500 = 1.0f/500.0f;
    unsigned int original, counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    char addRemoveText[50];

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

//    double tempContactSum, currentContactsSum = contactPotential(&beads_in_use_tree, pModel);
//    double tempKLDivContacts, currentKLDivContacts = currentContactsSum/(double)beads_in_use_tree.size();
    totalViolations = getViolationsFromSet(&beads_in_use_tree, contactsDistributionOfModel, workingLimit, pModel);
//    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    double tempKLDivContacts = 0, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyFivePercent = (unsigned int)(0.65*step_limit);
    const auto eightyfivePercent = (unsigned int)(0.85*step_limit);
    const auto ninetyfivePercent = (unsigned int)(step_limit - 3*workingLimit);

    const auto totalSubunits = (double)pModel->getNumberOfSubUnits();

    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = (isRefine) ? 0 : currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
    double contactsAvg = 2.0*binCount[0]/((float)workingLimit*totalSubunits);

    double this_energy, current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
    // Add violations
    current_energy += beta*totalViolations;

    double lowest_energy = current_energy ;

    logger("alpha", formatNumber(alpha, 3));
    logger("CDW(alpha)", formatNumber(contactDistributionWeight, 3));
    logger("lambda", formatNumber(lambda, 3));
    logger("beta", formatNumber(beta, 3));

    logger("NUMBER OF COMP", std::to_string(currentNumberOfComponents));
    logger("VIOLATIONS", std::to_string(totalViolations));
    logger("D_KL", formatNumber(currentKL, 7));
    logger("DIV_CON", formatNumber(currentKLDivContacts, 4));
    logger("TOTAL ENERGY", formatNumber(current_energy, 5));

    std::clock_t startTime;
    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    std::vector<unsigned int> selections(workingLimit);
    for(unsigned int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }


    unsigned int numberOfCoolingTempSteps=0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        auto set_it = hull.begin();
        auto beginBinCount = binCount.begin();
        startTime = std::clock();

        if ( distribution(gen) < percentAddRemove ){ //add or remove bead within working Set (exclude deadzone)
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            if (distribution(gen) < 0.51){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                // find a bead to add to that is weighted by number of contacts it has
                original = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                while(original == pModel->getNeighborLimit()){
                    original = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                }

//                std::advance(set_it, randomHull(gen));
//                original  = *set_it;
                auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

                std::iter_swap(bead_indices.begin() + workingLimit, itIndex); // this swaps to working position and changes backup
                // increment workingLimit to include new position
                workingLimit += 1;

                beads_in_use_tree.insert(original);
                addToPrSym(original, bead_indices, workingLimit, binCount, pModel, pData);
                testKL = pData->getScore(binCount);

                tempNumberOfComponents = eulerTour.addNode(original, pModel);

                temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

//                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

                if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     ADD => %i", 1);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased
                } else { // undo changes (rejecting)
                    eulerTour.removeNode(original);
                    beads_in_use_tree.erase(original);
                    std::sprintf(addRemoveText, "     ADD => %i", 0);
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?

                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                original = bead_indices[selections[0]];
                // test for deletion
                bool foundIt = true;
                if (numberOfCoolingTempSteps > ninetyfivePercent){
                    for(const auto & select : selections){
                        original = bead_indices[select];
                        if (numberOfContactsFromSet(&beads_in_use_tree, pModel, original) == 1){
                            tempNumberOfComponents = eulerTour.removeNode(original);
                            foundIt = false;
                            break;
                        }
                    }
                } else {
                    for(const auto & select : selections){
                        tempNumberOfComponents = eulerTour.removeNode(original);
                        if (tempNumberOfComponents <= currentNumberOfComponents){
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(original, pModel);
                        original = bead_indices[select];
                    }
                }

                if (foundIt){
                    original = bead_indices[randomIndex(gen)];
                    tempNumberOfComponents = eulerTour.removeNode(original);
                }

                removeLatticePositionToModelSym(bead_indices, binCount, &workingLimit, &original, pModel, pData);
                beads_in_use_tree.erase(original); // remove bead for next calculation

                temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
                testKL = pData->getScore(binCount);

                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

//                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

                if (this_energy < current_energy || (std::exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     REM => %i", original);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased
                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
//                    auto pBeginIt = bead_indices.begin();
//                    restoreRemovingLatticePointFromBackUp(&pBeginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    workingLimit += 1;
                    // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
                    // sorting is n*log(n) with n = 200?  should be much smaller
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                    eulerTour.addNode(original, pModel);
                }
            }

        } else { // positional refinement
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            bool retry = true;
            if (numberOfCoolingTempSteps > ninetyfivePercent){
                for(const auto & select : selections){
                    position = select;
                    swap1 = bead_indices[position];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1) < 2){
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            retry = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if (retry){
                bool foundIt = true;
                for(const auto & select : selections){ // select bead to move
                    position = select;
                    swap1 = bead_indices[position];
                    if ((1-getProbabilityOfSelectedIndex(&beads_in_use_tree, pModel, swap1)) < distribution(gen)){ // favors beads with greater contacts
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            retry = false;
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }

                if(foundIt){
                    position = selections[randomIndex(gen)];
                    swap1 = bead_indices[position];
                    eulerTour.removeNode(swap1);
                }
            }
//
//
//            if (retry){ // favors densly packed areas
//                for(const auto & select : selections){
//                    position=select;
//                    swap1 = bead_indices[position];
//                    if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
//                        foundIt = "true";
//                        break;
//                    }
//                    eulerTour.addNode(swap1, pModel);
//                }
//            }

            //std::cout << "selected " << foundIt << std::endl;
            // remove contribution of swap1
            auto itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
            beads_in_use_tree.erase(swap1);

            // find better position weighted by the number of contacts
//            std::advance(set_it, randomHull(gen));
//            unsigned int originalSwap2Value = *set_it;
//            while (originalSwap2Value == swap1){
//                set_it = hull.begin();
//                std::advance(set_it, randomHull(gen));
//                originalSwap2Value = *set_it;
//            }

            unsigned int originalSwap2Value = pModel->getNeighborLimit(); // find new position weighted by connectivity
            for(auto & select : selections){ // array indices max value is workingLimit
                originalSwap2Value = getUseableNeighborWeighted(&beads_in_use_tree, pModel, bead_indices[select]);
                if (originalSwap2Value != pModel->getNeighborLimit() && originalSwap2Value != swap1){
                    break;
                }
            }

            // get available neighbor position
            auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);
            //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted
            addToPrSym(originalSwap2Value, bead_indices, workingLimit, binCount, pModel, pData);

            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);
            beads_in_use_tree.insert(originalSwap2Value);

//            tempContactSum = addToContactsPotential(originalSwap2Value, tempContactSum, &beads_in_use_tree, pModel);
//            tempKLDivContacts=tempContactSum/(double)beads_in_use_tree.size();
            temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
//            populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
//            subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//            temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;
            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

            if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)) ) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
            } else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                beads_in_use_tree.erase(originalSwap2Value);
                std::sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        printf("      TEMP :  %-.3E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %10.5f FAILURES => %6i  TIME : %.5f\n", acceptRate, failures, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
        printf("     LIMIT : %10i     COMP => %3i VIOL => %d EVOL %.2E \n", workingLimit, currentNumberOfComponents, totalViolations, current_volume_energy);
        printf("  CONTACTS :   %.2E (%.2E ALPHA : %.1E) <AVE> : %.2f  \n", (currentKLDivContacts*contactDistributionWeight), currentKLDivContacts, currentCDW, contactsAvg);
        printf("       SYM : %10s  OUTFILE => %s\n", sym.c_str(), pnameTo);
        printf("      D_KL :  %-5.3E ( %.3E ) TOTAL ENRGY : %.4E\n", currentKL, lowest_energy, current_energy);

        // update run time parameters
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = totalViolations;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            currentKLDivContacts = tempKLDivContacts;
            current_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

//            if (numberOfCoolingTempSteps > sixtyFivePercent){
//            populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
//            } else {
//                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
//            }
//            randomHull = std::uniform_int_distribution<unsigned int>(0,(unsigned int)hull.size()-1);
//            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/((float)workingLimit*totalSubunits);

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;

//        calculateModelPrDistributionSym(&bead_indices, &testBinCount, workingLimit, temp_violations, pModel, pData );
//        float testKL1 = pData->getScore(binCount);
//        if ((currentKL != testKL1)  || checkForRepeats(bead_indices)){
//            testKL = pData->getScore(binCount);
//            std::cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << " -- " << testKL << std::endl;
//            exit(0);
//        }

        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent ){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy; // + totalContactEnergy ;
        }

        if (once && numberOfCoolingTempSteps > eightyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy;
            once = false;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        counter++;
    } // end of steps


    calculateModelParametersSymmetry(&beads_in_use_tree, pModel);

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
        unsigned int temp = bead_indices[i];
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, temp);
    }

    tempAverageContacts = (float)(tempAverageContacts/(double)workingLimit);
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    // calculate Rg
    //std::cout << " CALCULATED RG : " << calculateRgSym(bead_indices, workingLimit, pModel) << std::endl;

    if (currentNumberOfComponents > 1 ) {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return "failed";
    }

    /*
     * calculate distribution of contacts
     */
    unsigned int totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    logger("PRINTING", "DISTRIBUTION OF CONTACTS");
    for (unsigned int c=1; c<13; c++){
        unsigned int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }

//    float volume_test = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//    pData->printKLDivergence(binCount);
    return nameOfModel;
}


/**
 * Model is remapped to existing lattice defined by bin_width in Data
 *
 * @param pModel
 * @param pData
 * @param name
 * @param PDBFilename, should be something like damstart/damfilt
 * @return
 */
bool Anneal::initializeModelToRefineSym(PointSetModel *pModel, PofRData *pData, std::string name, std::string PDBFilename) {

    srand(time(0));
    contactCutOff = pModel->getNeighborCutOffLimit();
    violation_limit = pModel->getBeadRadius();
    // float * pDistance = pModel->getPointerToDistance();
    // convert distances within the large search space to bins based on input P(R)-DATA file
    // this->fillPrBinsAndAssignTotalBin( pModel, pData);
    maxbin = pModel->getMaxBin(pData);
    // adjust this part
    maxbin += 2;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // initialize Universe and fill indices
    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's

    logger("TOTAL PTS IN UNIVERSE", std::to_string(totalBeadsInSphere));
    unsigned int * ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    // returns only seed, needs to be mapped to entire Universe
    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDBSym(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
    unsigned int workingLimit = pModel->getTotalInSeed();

    // sort reassigned bead model into working universe
    unsigned int locale=0;
    for(auto  it = pModel->getSeedBegin(); it != pModel->getSeedEnd(); ++it) {
        auto fit = std::find(bead_indices.begin(), bead_indices.end(), *it); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(bead_indices.begin()+locale, fit);
        locale++;
    }

    if (locale != workingLimit){
        std::cout << " Locale does not equal working limit " << std::endl;
        exit(0);
    }

    std::sort(bead_indices.begin(), bead_indices.begin() + locale);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50

    logger("EXP N_S BINS", std::to_string(totalBins));
    logger("MODEL N_S BINS", std::to_string(maxbin));
    logger("BINWIDTH", formatNumber((float)pData->getBinWidth(), 2));
    logger("BEAD RADIUS", formatNumber(pModel->getBeadRadius(), 2));

    unsigned int violations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, violations, pModel, pData );
    // fill binCount for first time
    auto currentKL = (float)pData->getScore(binCount);
    //
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    logger("WORKINGLIMIT", std::to_string(workingLimit));
    logger("VIOLATIONS", std::to_string(violations));
    // parse violations?
    // setup parameters for hull
    // area of the model P(r) distribution
    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents  = eulerTour.getNumberOfComponents();

    if (currentNumberOfComponents > 1){
        logger("ERROR UNFIT STARTING MODEL", PDBFilename);
        return false;
    }
    logger("STARTING D_KL", formatNumber(currentKL, 7));
    char flags[] = "qhull FA";
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // int high;
    // pModel->writeModelToFile(workingLimit, bead_indices, name, high);
    pModel->setReducedSeed(workingLimit, bead_indices);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    pModel->setBeadAverageAndStdev(workingLimit, 0);

    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                 VOLUME                 *******************" << std::endl;
    logger("CONVEX HULL VOLUME", formatNumber(current_volume,1));
    std::cout << "*******************                                        *******************" << std::endl;

    logger("SHANNON BINS IN DATA", std::to_string(pData->getShannonBins()));
    logger("ZERO BIN", std::to_string( pData->getZeroBin()));
    logger("","EXITING INITIAL MODEL BUILD");

    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        unsigned int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }

//    double runningContactsSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    logger("AVERAGE NUMBER CONTACTS", formatNumber(average_number_of_contacts,1));
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_remapped_x", this, pData, 0, current_volume, average_number_of_contacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_remapped_x", this, pData, 0, current_volume, average_number_of_contacts);

    isRefine = true;
    return true;
}


/**
 * Input model must be aligned to a subunit
 * Input model is the masked Universe
 * Starting from a filled space, sample down
 * Starting Model has a connected Euler tour
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineSymModelRefine(PointSetModel *pModel, PofRData *pData, std::string nameTo){
    const std::string sym = pModel->getSymmetry();
    bool once = true;
    const char * pnameTo = nameTo.c_str() ;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.0000001 : (double)highTempStartForCooling;
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;

    //unsigned int workingLimit = pModel->getStartingWorkingLimit(); // restrict universe to size of the input model that is going to be refined
    //const unsigned int totalBeadsInSphere = workingLimit;
    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    /*
     * copy model into bead_indices
     */
    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit(); // restrict universe to size of the input model that is going to be refined

    // randomize upto workingLimit and use original
    std::shuffle(bead_indices.begin(), bead_indices.begin()+workingLimit, gen);
    workingLimit = pModel->getBaseWorkingLimit();
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    //workingLimit = minimizeViolations(&bead_indices, workingLimit, pModel);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin()+workingLimit);

    /*
     * create set of beads to select from based on all those that are in contact with selected set
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);

    if (maxbin < 2){
        throw std::invalid_argument(" INCORRECT MAXBIN NOT INITIALIZED PROPERLY");
    }

    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin);   // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);      // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", std::to_string(totalBins));
    logger("MAX MODEL N_S BINS", std::to_string(totalBins));
    logger("BINWIDTH (Angstroms)", formatNumber((float)pData->getBinWidth(), 2));

    unsigned int temp_violations, totalViolations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, totalViolations, pModel, pData );
    // fill binCount for first time
    double testKL, currentKL = pData->getScore(binCount);
    char flags[] = "qhull FA"; // CVX HULL STUFF

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH" << std::endl;
    logger("STARTING D_KL", formatNumber(currentKL, 5));

    float tempAverageContacts;
    // coupon collector's problem
    const unsigned int updateCount = ccmultiple*getCoupons(workingLimit);

    float step_limit = (updateCount < 51357) ? 51357 : (float)updateCount;

    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int tempNumberOfComponents = 0, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool readyToAdd = true;
//    auto readyToAddLimit = (unsigned int)(totalBeadsInSphere*0.71);
    auto readyToAddLimit = pModel->getBaseWorkingLimit();


    bool isUpdated = false;
    float acceptRate = 0.5;
    const float inv500slash499 = 499.0f/500.0f, inv500 = 1.0f/500.0f;
    unsigned int original, counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    char addRemoveText[50];

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

//    double tempContactSum, currentContactsSum = contactPotential(&beads_in_use_tree, pModel);
//    double tempKLDivContacts, currentKLDivContacts = currentContactsSum/(double)beads_in_use_tree.size();
    totalViolations = getViolationsFromSet(&beads_in_use_tree, contactsDistributionOfModel, workingLimit, pModel);
    double tempKLDivContacts = 0, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    auto fifteenPercent = (unsigned int)(0.15*step_limit);
    auto sixtyFivePercent = (unsigned int)(0.65*step_limit);
    auto eightyfivePercent = (unsigned int)(0.85*step_limit);

    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = (isRefine) ? 0 : currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
    double contactsAvg=3.0d;

    double this_energy, current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;

    double lowest_energy = current_energy ;

    logger("alpha", formatNumber(alpha, 3));
    logger("CDW(alpha)", formatNumber(contactDistributionWeight, 3));
    logger("lambda", formatNumber(lambda, 3));
    logger("beta", formatNumber(beta, 3));

    logger("NUMBER OF COMP", std::to_string(currentNumberOfComponents));
    logger("VIOLATIONS", std::to_string(totalViolations));
    logger("D_KL", formatNumber(currentKL, 7));
    logger("DIV_CON", formatNumber(currentKLDivContacts, 4));
    logger("TOTAL ENERGY", formatNumber(current_energy, 5));

    std::clock_t startTime;
    unsigned int hullsize = hull.size();
    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    std::vector<unsigned int> selections(workingLimit);
    unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
    for(unsigned int i=0; i < workingLimit; i++){
        *(pSelections+i) = i; // some distances will exceed dmax
    }

    const auto totalSubunits = (double)pModel->getNumberOfSubUnits();
    unsigned int numberOfCoolingTempSteps=0;
    bool doIt=false;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        startTime = std::clock();

        std::cout << "______________________________________________________________________________" << std::endl;
        std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
        // additional points to expand deadlimit will occur via enlarging CVX Hull
        if (readyToAdd && distribution(gen) < 0.499){ // ADD BEAD?
            std::cout << "*******************                  ADD                   *******************" << std::endl;
            auto set_it = hull.begin();
            std::advance(set_it, randomHull(gen));
            original  = *set_it;
            auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

            std::iter_swap(bead_indices.begin() + workingLimit, itIndex); // this swaps to working position and changes backup
            // increment workingLimit to include new position
            workingLimit += 1;

            beads_in_use_tree.insert(original);
            addToPrSym(original, bead_indices, workingLimit, binCount, pModel, pData);
            testKL = pData->getScore(binCount);

            tempNumberOfComponents = eulerTour.addNode(original, pModel);

            temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

            if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     ADD => %i", 1);
                hull.erase(original);
            } else { // undo changes (rejecting)
                eulerTour.removeNode(original);
                beads_in_use_tree.erase(original);
                std::sprintf(addRemoveText, "     ADD => %i", 0);
                auto beginBinCount = binCount.begin();
                restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
            }

        } else { // REMOVE BEADS?

            std::cout << "*******************                 REMOVE                 *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
            //original = bead_indices[randomIndex(gen)];
            original = bead_indices[selections[0]];
            // test for deletion
            if (numberOfCoolingTempSteps > sixtyFivePercent){
                for(const auto & select : selections){
                    original = bead_indices[select];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, original) < contactsAvg){

                        tempNumberOfComponents = eulerTour.removeNode(original);
                        if (tempNumberOfComponents <= currentNumberOfComponents){
                            break;
                        }
                        eulerTour.addNode(original, pModel);
                    }
                }
            } else {
                for(const auto & select : selections){
                    tempNumberOfComponents = eulerTour.removeNode(original);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[select];
                }
            }

            removeLatticePositionToModelSym(bead_indices, binCount, &workingLimit, &original, pModel, pData);
            beads_in_use_tree.erase(original); // remove bead for next calculation

            temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
            testKL = pData->getScore(binCount);

            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

            if (this_energy < current_energy || (std::exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) )) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     REM => %i", original);
                hull.insert(original);
            } else { // undo changes and move to next bead (rejecting)
                beads_in_use_tree.insert(original);
//                    auto pBeginIt = bead_indices.begin();
//                    restoreRemovingLatticePointFromBackUp(&pBeginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                workingLimit += 1;
                // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
                // sorting is n*log(n) with n = 200?  should be much smaller
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                eulerTour.addNode(original, pModel);
            }
        }

        printf("      TEMP :  %-.3E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %10.5f FAILURES => %6i  TIME : %.5f\n", acceptRate, failures, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
        printf("       SYM : %10s  OUTFILE => %s\n", sym.c_str(), pnameTo);
        printf("    E_VIOL :   %.2E VIOL => %d LIMIT : %5i COMP %4i HULL %i\n", (beta*totalViolations), totalViolations, workingLimit, currentNumberOfComponents, hullsize);
        printf("  E_CNTCTS :   %.2E (%.2E ALPHA : %.1E) <AVE> : %.2f  \n", (currentKLDivContacts*contactDistributionWeight), currentKLDivContacts, currentCDW, contactsAvg);
        printf("      D_KL :  %-5.3E ( %.3E ) TOTAL ENRGY : %.4E\n", currentKL, lowest_energy, current_energy);

        // update run time parameters
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = totalViolations;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            currentKLDivContacts = tempKLDivContacts;
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            /*
             * do not update hull until target acceptance rate is reached
             */
            if(doIt){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
            }
            hullsize = hull.size();
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);
//            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1);
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

            if (workingLimit <= readyToAddLimit){
                readyToAdd=true;
            }
//            currentContactsSum = tempContactSum;
            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/((float)workingLimit*totalSubunits);
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        if (!doIt && numberOfCoolingTempSteps > fifteenPercent){
            doIt = true;
        }

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;
//        calculateModelPrDistributionSym(&bead_indices, &testBinCount, workingLimit, temp_violations, pModel, pData );
//        double testKL1 = pData->getScore(binCount);
//        if ((currentKL != testKL1)  || checkForRepeats(bead_indices)){
//            testKL = pData->getScore(binCount);
//            std::cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << " -- " << testKL << std::endl;
//            exit(0);
//        }

        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent ){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
//            current_energy = currentKL + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts;
            lowest_energy = current_energy;// + totalContactEnergy ;
        }

        if (once && numberOfCoolingTempSteps > eightyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts;
//            current_energy = currentKL + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy;
            once = false;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        counter++;
    } // end of steps


    calculateModelParametersSymmetry(&beads_in_use_tree, pModel);

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
        unsigned int temp = bead_indices[i];
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, temp);
    }

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    tempAverageContacts = (float)(tempAverageContacts/(double)workingLimit);
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    if (currentNumberOfComponents > 1 ) {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return "failed";
    }

    /*
     * calculate distribution of contacts
     */
    unsigned int totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (unsigned int c=1; c<13; c++){
        unsigned int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }


    // At end of each temp, update a probability model for volume?  Use this to select
    // perform positional refinement until delta E stabilizes?
    // final round to move any points close to body
    float volume_test = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    // pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
    // pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");

//    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_final_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);
//    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);

    //pModel->writeModelToFile(workingLimit, bead_indices, "refined_");
    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");

    //  pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_" + filenameprefix);
    //  pModel->writeSymModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_sym_");

    return nameOfModel;
}