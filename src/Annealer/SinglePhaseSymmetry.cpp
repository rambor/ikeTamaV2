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
 * Speed changes (all exact):
 *  - symmetry mates of every lattice point are computed once (ensureSymmetryCache); P(r) updates use them
 *  - clashes with symmetry mates are counted incrementally (ClashTracker) instead of O(S*N^2) per move
 *  - ContactTracker for occupancy / contact numbers / contact energy, SwapSelection for bead_indices
 *  - hull points drawn from a vector snapshot instead of std::advance over a std::set for every draw
 *  - convex hull volume recomputed only when the hull can change (HullVolumeCache)
 *  - progress printed every few hundred steps instead of every step
 *
 * refineSymModel / refineSymModelRefine also contain bug fixes, listed at each function.
 */

#ifdef ANNEAL_VERIFY
// independent clash count: every bead against every symmetry mate (the original getViolations)
static unsigned int bruteForceViolations(const std::vector<unsigned int> & v, unsigned int wl, PointSetModel * pModel, float limit){
    const unsigned int S = pModel->getNumberOfSubUnits();
    std::vector<vector3> coordinates(S*wl);
    for (unsigned int i=0; i<wl; i++) coordinates[i] = pModel->getBead(v[i])->getVec();
    unsigned int count = wl;
    for (unsigned int s=1; s<S; s++) pModel->transformCoordinatesBySymmetryPreCalc(s, wl, count, coordinates);
    unsigned int violations = 0;
    for (unsigned int i=0; i<wl; i++)
        for (unsigned int t=wl; t<S*wl; t++)
            if ((coordinates[i] - coordinates[t]).length() < limit) violations++;
    return violations;
}
#endif


/**
 * create initial model of complete object
 * does not require pre-computed distances, this is a direct method meaning all pairwise distances are calculated on the fly
 */
bool Anneal::createInitialModelSymmetry(PointSetModel *pModel, PofRData *pData) {

    srand(time(0));
    std::mt19937 gen(seedForRun(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    auto lowTempStop =  (double)0.0001;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    contactCutOff = pModel->getNeighborCutOffLimit();
    violation_limit = pModel->getBeadRadius();

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
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> tempContactsDistributionOfModel(13);

    logger("TOTAL EXP N_S BINS", formatNumber(totalBins));
    logger("MAX MODEL N_S BINS", formatNumber(maxbin));
    logger("BINWIDTH", formatNumber((float)pData->getBinWidth(), 2));
    logger("BEAD RADIUS", formatNumber((float)pModel->getBeadRadius(), 2));

    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // as bead indices are discarded, set upper limit of vector
    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> lowest_subUnit_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's

    const unsigned int num = totalBeadsInSphere;
    unsigned int * const ptr = (num != 0) ? subUnit_indices.data() : nullptr;
    for(unsigned int i = 0; i < num; i++) {
        ptr[i] = i;
    }

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.331);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    const double targetSubUnitVolume = (2*lowerV+upperV)/3/(float)pModel->getNumberOfSubUnits();
    const auto radius_larger = (float)(pModel->getBeadRadius()*std::sqrt(7.0)/2.0);
    const auto lowerN = (unsigned int)(std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/3;
    const auto upperN = (unsigned int)(std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()));

    // the bead count comes from the Porod volume in the P(r) file; a missing or zero VOLUME gives no beads, and qhull
    // then fails on an empty point set (segmentation fault)
    if (pData->getVolume() <= 0 || (upperN + lowerN)/2 < 4 || upperN > totalBeadsInSphere){
        std::ostringstream msg;
        msg << "\n  ERROR createInitialModelSymmetry: cannot size the initial model from the P(r) file " << pData->getFilename()
                  << "\n        VOLUME : " << pData->getVolume() << " A^3 gives " << (upperN + lowerN)/2 << " beads per subunit of radius "
                  << pModel->getBeadRadius() << " (need 4 to " << totalBeadsInSphere << ").\n        Set a Porod volume (VOLUME) in the P(r) file.";
        throw std::invalid_argument(msg.str());
    }
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
    std::vector<unsigned int> hullPoints(hull.begin(), hull.end());   // same order as the set, O(1) random access
    std::string switched = "HULL";
    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    std::vector<coordT> points(3*(size_t)totalBeadsInSphere);
    char flags[25];
    std::sprintf(flags, "qhull s FA");

    HullVolumeCache hullVolume(pModel, flags, subUnit_indices, subUnitWorkingLimit);
    float subunit_test_volume, subunit_volume = hullVolume.currentVolume();

    // occupancy of the accepted set (mirrors beads_in_use_tree) and bookkeeping
    ContactTracker occupancy(pModel->getDirectPointerToNeighborhood(), pModel->getSizeOfNeighborhood(), pModel->getNeighborLimit(), totalBeadsInSphere);
    occupancy.reset(subUnit_indices.begin(), subUnitWorkingLimit);
    SwapSelection selection(subUnit_indices, subUnitWorkingLimit, totalBeadsInSphere);
    ensureSymmetryCache(pModel);
    ensureClashTable(pModel);
    ifaceHard = false;   // initial model: soft requirement with a gap gradient,
    ifaceArmed = false;  // switched on once the subunit is first one connected piece (does not compete with connectivity)
    trkReset(pModel, subUnit_indices.begin(), subUnitWorkingLimit);   // trial set

    // calculate starting energy
    unsigned int totalViolations, temp_violations;
    std::clock_t start = std::clock();
    calculateModelPrDistributionSym(&subUnit_indices, &binCount, subUnitWorkingLimit, totalViolations, pModel, pData );

    // fill binCount for first time
    double currentKL = pData->getScore(binCount);

    double hlambda = pData->getIsPr() ? 0.01 : 10.0;
    double muPercent = 0.731;
    double muConstant = currentKL*muPercent/((1.0 - muPercent)*subunit_volume/(double)subUnitWorkingLimit);

    unsigned int failures=0;

    bool updateCVX = false;

    logger("THREAD TIME (sec)", formatNumber( (float)((std::clock() - start)/(double) CLOCKS_PER_SEC), 6) );

    double testKL, test_energy; // sets alpha as a constant during high temp eq

    unsigned int lowestWorkingLimit = subUnitWorkingLimit;
    unsigned int tempNumberOfComponents, currentNumberOfComponents;

    EulerTour eulerTour(subUnit_indices.begin(), subUnitWorkingLimit, pModel);
    currentNumberOfComponents = eulerTour.getNumberOfComponents();

    double current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) ;
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
    current_energy += beta*totalViolations + ifacePenalty();

    double lowest_energy = current_energy;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(binCount.begin(), binCount.end(), lowestbinCount.begin());  // unaltered P(r)

    unsigned int high, hullsize = (unsigned int)hull.size();
    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased

    double temp_subunit_volume_energy, current_subunit_volume_energy = muConstant*subunit_volume/(double)subUnitWorkingLimit;

    unsigned int addFailed = 0, removeFailed =0;
    const unsigned int printEvery = 397;
    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        if (distribution(gen) < 0.17 || currentNumberOfComponents > 5){
            // create points from workSet to determine HULL
            std::sprintf(addRemoveText, "POSITIONAL");
            unsigned int position = randomIndex(gen);
            unsigned int swap1 = subUnit_indices[position];

            // setup parameters for hull
            if ( (distribution(gen) < 0.37129) || currentNumberOfComponents > 2 ){ // only move CVX points if true

                unsigned int * const pSelections = subUnit_indices.data(); // initialized as empty in Model class
                for (unsigned int i = 0; i < subUnitWorkingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                    active_indices[i] = pSelections[i];
                }

                qh_new_qhull(3, subUnitWorkingLimit, points.data(), 0, flags, nullptr, nullptr);
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
                position = selection.position(swap1);
            } // recalculate

            // find bead to swap in active set
            removeFromPrSymAt(position, subUnit_indices, subUnitWorkingLimit, binCount, pData); // remove selected index from P(r)
            eulerTour.removeNode(swap1);
            // find position to swap to
            unsigned int neighbor = hullPoints[randomHull(gen)];
            while (neighbor == swap1 || occupancy.occupied(neighbor)) { // find a new position
                neighbor = hullPoints[randomHull(gen)];
            }

            // make the swap and update P(r)
            selection.swapAt(position, neighbor);
            addToPrSymAt(position, subUnit_indices, subUnitWorkingLimit, binCount, pData);
            trkRemove(swap1); trkAdd(neighbor);
            temp_violations = clashes.count();

            testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            subunit_test_volume = hullVolume.afterSwap(swap1, neighbor, subUnit_indices, subUnitWorkingLimit); // calculate volume of the subUnit
            temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;

            test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) +  beta*temp_violations + ifacePenalty();//diff*invTarVolmuConstant;

            if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                beads_in_use_tree.erase(swap1);
                beads_in_use_tree.insert(neighbor);
                occupancy.remove(swap1); occupancy.add(neighbor);
                hullVolume.accept();
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                selection.undoSwap();
                trkRemove(neighbor); trkAdd(swap1);
                hullVolume.reject();
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

                // == removeLatticePositionToModelSym
                removeFromPrSymAt(selection.position(testIndex), subUnit_indices, subUnitWorkingLimit, binCount, pData);
                selection.remove(testIndex);
                trkRemove(testIndex);
                temp_violations = clashes.count();

                testKL = pData->getScore(binCount);

                subunit_test_volume = hullVolume.afterRemove(testIndex, subUnit_indices, subUnitWorkingLimit); // calculate volume of the subUnit
                temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + ifacePenalty();

                if ( ( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy - temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.erase(testIndex);
                    occupancy.remove(testIndex);
                    hullVolume.accept();
                    updateCVX=true;
                } else { // undo changes and move to next bead (rejecting)
                    eulerTour.addNode(testIndex, pModel);
                    selection.undoRemove();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); //copy to bin count
                    trkAdd(testIndex);
                    hullVolume.reject();
                    removeFailed += 1;
                }

            } else { // ADD beads
                std::sprintf(addRemoveText, "    ADD   ");

                // if we select a random bead to add to, we will bias to, comparatively, denser regions
                unsigned int testIndex = getUseableNeighborWeighted(occupancy, pModel, subUnit_indices[randomIndex(gen)], gen);
                while(testIndex == pModel->getNeighborLimit()){
                    testIndex = getUseableNeighborWeighted(occupancy, pModel, subUnit_indices[randomIndex(gen)], gen);
                }

                selection.add(testIndex);   // swaps to working position, increments workingLimit

                addToPrSymAt(subUnitWorkingLimit-1, subUnit_indices, subUnitWorkingLimit, binCount, pData);
                trkAdd(testIndex);
                temp_violations = clashes.count();

                testKL = pData->getScore(binCount);

                subunit_test_volume = hullVolume.afterAdd(testIndex, subUnit_indices, subUnitWorkingLimit); // calculate volume of the subUnit
                temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel); // want a single tour (so implies connectivity)

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations + ifacePenalty();// + diff*invTarVolmuConstant + beta*temp_violations + ifacePenalty();

                if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.insert(testIndex);
                    occupancy.add(testIndex);
                    hullVolume.accept();
                    updateCVX=true;
                } else { // undo changes (rejecting)
                    selection.undoAdd();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    trkRemove(testIndex);
                    hullVolume.reject();
                    eulerTour.removeNode(testIndex);
                    addFailed += 1;
                }
            }
        }


        if (updateCVX){

            current_energy = test_energy;
            currentKL = testKL;
            currentNumberOfComponents = tempNumberOfComponents;
            if (ifaceOn && !ifaceArmed && currentNumberOfComponents == 1){ ifaceArmed = true; current_energy += ifacePenalty(); test_energy = current_energy; lowest_energy = current_energy; }
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
            hullPoints.assign(hull.begin(), hull.end());
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);

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

        if (high % printEvery == 0){
            std::printf("*******************             %s                 ******************* \n", addRemoveText);
            std::printf("      MAXSTEPS : %7i (%7i) ACCEPTRATE : %.3f TEMP : %.2E\n", highTempRounds, high, acceptRate, lowTempStop);
            std::printf("         GRAPH : %7i  MODE => %s \n", currentNumberOfComponents, switched.c_str());
            std::printf("         LIMIT : %7i UPPER <= %i  LOWER >= %i (%i)\n", subUnitWorkingLimit, upperN, lowerN, hullsize);
            std::printf("        VOLUME : %.0f (%.0f)  MU => %.1E  MU*VOL => %.2E\n", subunit_volume, targetSubUnitVolume, muConstant, current_subunit_volume_energy);
            std::printf("    VIOLATIONS : %7d  BETA => %.1E (%i) %i | %i  INTERFACE %u (min %u)\n", totalViolations, beta, volumeCount, addFailed, removeFailed, interfaceContacts(), ifaceMin);
            std::printf("          D_KL : %.4E ENRGY: %.4E (%.4E) \n", currentKL, current_energy, lowest_energy);
            std::cout << "*******************                                        *******************" << std::endl;
        }

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
            current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + beta*totalViolations + ifacePenalty();
        }

#ifdef ANNEAL_VERIFY
        {
            std::vector<unsigned int> fullBins(maxbin); unsigned int v;
            calculateModelPrDistributionSym(&subUnit_indices, &fullBins, subUnitWorkingLimit, v, pModel, pData);
            if (fullBins != binCount || bruteForceViolations(subUnit_indices, subUnitWorkingLimit, pModel, violation_limit) != totalViolations || eulerTour.validateList("initsym")){
                std::cout << "VERIFY FAILED (createInitialModelSymmetry) at step " << high << std::endl; exit(1);
            }
        }
#endif
    } // end of HIGH TEMP EQUILIBRATION

    logger("HULL VOLUMES COMPUTED", formatNumber(hullVolume.recomputed));
    logger("HULL VOLUMES REUSED", formatNumber(hullVolume.reused));

    switched = "HULLFinal";
    recalculateDeadLimit(lowestWorkingLimit, lowest_subUnit_indices, hull, pModel);
    pModel->writeSetToFile(hull, switched);


    highTempStartForCooling = (float)lowTempStop;
    // calculate average volume and standard deviation
    // no connected state was ever recorded (volumeCount = 0): use the current size instead of dividing by zero
    float volumeAverage = (volumeCount > 0) ? workingLimitSum/(float)volumeCount : (float)subUnitWorkingLimit;
    float volumeStdev = (volumeCount > 0) ? std::sqrt(std::max(0.0f, sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage)) : 0.0f;
    if (volumeCount == 0){
        std::cout << "  WARNING: the subunit never became one connected piece during the high-temperature search;" << std::endl;
        std::cout << "           increase highTempRounds (C6 with a 270 A Dmax needed about 10000)." << std::endl;
    }

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
 * Bug fixes relative to the previous version:
 *  1. after an accepted move, current_energy was overwritten with the contact energy
 *     ("current_energy = tempContactEnergy") and contactsEnergy was never updated
 *  2. in the positional move, "break" (when no bead could be picked by contact number) ended the whole annealing run
 *  3. an unused add position (addMe) was drawn on every ADD move
 *  4. the energy rebuilt when the contact weight is rescaled omitted beta*violations, so it was not comparable
 *     with the energies of the trial moves
 *  5. percentAddRemove (a member) was scaled down permanently during the run; a local copy is scaled now
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineSymModel(PointSetModel *pModel, PofRData *pData, std::string nameTo){
    const std::string sym = pModel->getSymmetry();
    const char * pnameTo = nameTo.c_str() ;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000000001 : (double)highTempStartForCooling;
    srand(time(0));

    std::string status = "RAMPING TO TEMP";

    std::mt19937 gen(seedForRun(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    const unsigned int noNeighborIndex = pModel->getNeighborLimit();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's

    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
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
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", std::to_string(totalBins));
    logger("MAX MODEL N_S BINS", std::to_string(totalBins));
    logger("BINWIDTH (Angstroms)", formatNumber(pData->getBinWidth(), 3));

    unsigned int temp_violations, totalViolations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, totalViolations, pModel, pData );
    // fill binCount for first time
    double testKL, currentKL = pData->getScore(binCount);

    char flags[] = "qhull FA"; // CVX HULL STUFF

    // based on the subunit estimated from initial modeling, use as a upper limit of the volume
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    const float target_volume = 0.8f*current_volume;
    double current_volume_energy = 0;

    logger("TARGET VOLUME (A^3)", formatNumber(target_volume, 3));
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

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

    // incremental bookkeeping
    ContactTracker contacts(pModel->getDirectPointerToNeighborhood(), pModel->getSizeOfNeighborhood(), noNeighborIndex, totalBeadsInSphere);
    contacts.reset(bead_indices.begin(), workingLimit);
    SwapSelection selection(bead_indices, workingLimit, totalBeadsInSphere);
    ensureSymmetryCache(pModel);
    ensureClashTable(pModel);
    trkReset(pModel, bead_indices.begin(), workingLimit);
    ifaceArmed = true; ifaceHard = ifaceOn && interfaces.count() >= ifaceMin;   // refinement: hard once the requirement is met
    totalViolations = clashes.count();

    // contacts potential
    double tempContactEnergy;
    double contactsEnergy = contacts.energy();

    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);
    const auto tenPercent = (unsigned int)(0.10*step_limit);
    const auto ninetyfivePercent = (unsigned int)(0.95*step_limit);

    float runPercentAddRemove = percentAddRemove;
    double alphaConstant = std::exp(1.0/(0.75*step_limit)*std::log(0.01/runPercentAddRemove));
    // scale percentAddRemove from 0.1 to 1 throughout the modeling

    double contactDistributionWeight = currentKL/contactsEnergy;
    // calculate and set starting energy of the model
    double current_energy = currentKL +
            lambda*connectivityPotential(currentNumberOfComponents) +
            contactDistributionWeight*contactsEnergy;

    double this_energy, lowest_energy = current_energy;

    // Add violations
    current_energy += beta*totalViolations + ifacePenalty();

    logger("alpha", formatNumber(alpha, 3));
    logger("lambda", formatNumber(lambda, 3));
    logger("beta", formatNumber(beta, 3));

    logger("NUMBER OF COMP", std::to_string(currentNumberOfComponents));
    logger("VIOLATIONS", std::to_string(totalViolations));
    logger("D_KL", formatNumber(currentKL, 7));
    logger("TOTAL ENERGY", formatNumber(current_energy, 5));

    std::clock_t startTime;
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    std::vector<unsigned int> selections(workingLimit);
    for(unsigned int i=0; i < workingLimit; i++){
        selections[i] = i;
    }

    auto rebuildEnergy = [&](){
        return currentKL + contactDistributionWeight*contactsEnergy
               + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + ifacePenalty();
    };

    unsigned int numberOfCoolingTempSteps=0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        startTime = std::clock();

        if ( distribution(gen) < runPercentAddRemove ){ //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            if (distribution(gen) < 0.4951){ // ADD BEAD?
                // find a bead to add to that is weighted by number of contacts it has
                original = getUseableNeighborWeighted(contacts, pModel, bead_indices[randomIndex(gen)], gen);
                while(original == noNeighborIndex){
                    original = getUseableNeighborWeighted(contacts, pModel, bead_indices[randomIndex(gen)], gen);
                }

                selection.add(original);                     // swaps to working position, workingLimit += 1
                contacts.add(original);
                addToPrSymAt(workingLimit-1, bead_indices, workingLimit, binCount, pData);
                testKL = pData->getScore(binCount);

                tempNumberOfComponents = eulerTour.addNode(original, pModel);

                trkAdd(original);
                temp_violations = clashes.count();
                tempContactEnergy = contacts.energy();

                this_energy = testKL +
                    lambda*connectivityPotential(tempNumberOfComponents) +
                        beta*temp_violations + ifacePenalty() +
                            contactDistributionWeight*tempContactEnergy;

                if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     ADD => %i", 1);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased
                } else { // undo changes (rejecting)
                    eulerTour.removeNode(original);
                    contacts.remove(original);
                    trkRemove(original);
                    std::sprintf(addRemoveText, "     ADD => %i", 0);
                    selection.undoAdd();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                }

            } else { // REMOVE BEADS?

                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                original = bead_indices[selections[0]];
                // test for deletion
                bool foundIt = true;
                if (numberOfCoolingTempSteps > ninetyfivePercent){
                    for(const auto & select : selections){
                        original = bead_indices[select];
                        if (contacts.contacts(original) == 1){
                            tempNumberOfComponents = eulerTour.removeNode(original);
                            foundIt = false;
                            break;
                        }
                    }
                } else {
                    for(const auto & select : selections){
                        original = bead_indices[select];
                        tempNumberOfComponents = eulerTour.removeNode(original);
                        if (tempNumberOfComponents <= currentNumberOfComponents){
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(original, pModel);
                    }
                }

                if (foundIt){
                    original = bead_indices[randomIndex(gen)];
                    tempNumberOfComponents = eulerTour.removeNode(original);
                }

                removeFromPrSymAt(selection.position(original), bead_indices, workingLimit, binCount, pData);
                selection.remove(original);                  // workingLimit -= 1
                contacts.remove(original); // remove bead for next calculation
                trkRemove(original);

                temp_violations = clashes.count();
                testKL = pData->getScore(binCount);

                tempContactEnergy = contacts.energy();

                this_energy = testKL +
                    lambda*connectivityPotential(tempNumberOfComponents) +
                        beta*temp_violations + ifacePenalty() +
                            contactDistributionWeight*tempContactEnergy;

                if (this_energy < current_energy || (std::exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     REM => %i", original);
                    randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased
                } else { // undo changes and move to next bead (rejecting)
                    contacts.add(original);
                    trkAdd(original);
                    selection.undoRemove();                  // workingLimit += 1
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    eulerTour.addNode(original, pModel);
                }
            }

        } else { // positional refinement
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            bool foundIt = true;
            if (numberOfCoolingTempSteps > sixtyfivePercent){
                for(const auto & select : selections){
                    position = select;
                    swap1 = bead_indices[position];
                    if (contacts.contacts(swap1) < 2){
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if (foundIt){
                for(const auto & select : selections){ // select bead to move
                    position = select;
                    swap1 = bead_indices[position];
                    auto tt = contacts.contacts(swap1);

                    if (tt < 5 || tt > 8) {
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            foundIt = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if (foundIt){ // nothing picked by contact number: first bead whose removal keeps the graph connected
                for(const auto & select : selections){
                    position = select;
                    swap1 = bead_indices[position];
                    if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                        foundIt = false;
                        break;
                    }
                    eulerTour.addNode(swap1, pModel);
                }
            }

            if (foundIt){ // every bead is a cut point; take the first one anyway
                position = selections[0];
                swap1 = bead_indices[position];
                eulerTour.removeNode(swap1);
            }

            // remove selected index from P(r)
            removeFromPrSymAt(position, bead_indices, workingLimit, binCount, pData);
            contacts.remove(swap1);
            trkRemove(swap1);

            unsigned int originalSwap2Value = noNeighborIndex; // find new position weighted by connectivity
            for(auto & select : selections){ // array indices max value is workingLimit
                originalSwap2Value = getUseableNeighborWeighted(contacts, pModel, bead_indices[select], gen);
                if (originalSwap2Value != noNeighborIndex && originalSwap2Value != swap1){
                    break;
                }
            }

            if (originalSwap2Value == noNeighborIndex || originalSwap2Value == swap1){ // nowhere to go: reject
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                contacts.add(swap1);
                trkAdd(swap1);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                std::sprintf(addRemoveText, "      FAILED => %i", swap1);
            } else {
                selection.swapAt(position, originalSwap2Value);
                addToPrSymAt(position, bead_indices, workingLimit, binCount, pData);

                // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
                testKL = pData->getScore(binCount);
                contacts.add(originalSwap2Value);
                trkAdd(originalSwap2Value);

                temp_violations = clashes.count();
                tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
                tempContactEnergy = contacts.energy();

                this_energy = testKL +
                    lambda*connectivityPotential(tempNumberOfComponents) +
                        beta*temp_violations + ifacePenalty() +
                            contactDistributionWeight*tempContactEnergy;

                if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)) ) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
                } else {
                    selection.undoSwap();
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                    contacts.remove(originalSwap2Value);
                    contacts.add(swap1);
                    trkRemove(originalSwap2Value);
                    trkAdd(swap1);
                    std::sprintf(addRemoveText, "      FAILED => %i", swap1);
                    eulerTour.removeNode(originalSwap2Value);
                    currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                }
            }
        } // end of positional refinement or add/remove if statement


        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            contactsEnergy = tempContactEnergy;            // (was: current_energy = tempContactEnergy)

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

        } else {
            if (!ifaceBlocked){ acceptRate = inv500slash499*acceptRate; }   // a move blocked by the interface constraint is not a failed trial
            if (!ifaceBlocked) failures++;
            ifaceBlocked = false;
        }

        // rescale etaFactor in small steps until final value
        // rescale contactsWeight in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent){
            status = "RAMPING";
            contactDistributionWeight = currentKL/contactsEnergy;//*0.733;
            current_energy = rebuildEnergy();
        } else if (numberOfCoolingTempSteps < sixtyfivePercent) {
            status = "EQUILIBRATING";
        }

        if (numberOfCoolingTempSteps < 0.75*step_limit) {
            runPercentAddRemove *= alphaConstant;
        }

        if (numberOfCoolingTempSteps > sixtyfivePercent && numberOfCoolingTempSteps%tenPercent == 0) { // reset the weight every 10%
            contactDistributionWeight = currentKL/contactsEnergy*0.5372;
            current_energy = rebuildEnergy();
        }

        if (numberOfCoolingTempSteps == sixtyfivePercent){
            contactDistributionWeight = currentKL/contactsEnergy*0.372;
            current_energy = rebuildEnergy();
        }

        if (current_energy < lowest_energy) lowest_energy = current_energy;

        if (numberOfCoolingTempSteps % 1001 == 0){
            printf("      TEMP : %-.3E MAXSTEPS => %.0f (%i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, status.c_str());
            printf("    ACCEPT : %10.5f FAILURES => %6i  TIME : %.5f\n", acceptRate, failures, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
            printf("     LIMIT : %10i     COMP => %3i VIOL => %d EVOL %.2E \n", workingLimit, currentNumberOfComponents, totalViolations, current_volume_energy);
            printf("  CONTACTS :   %.2E (%.2E )  \n", (contactsEnergy*contactDistributionWeight), contactsEnergy);
            if (ifaceOn && !ifaceHard && interfaces.count() >= ifaceMin) ifaceHard = true;
            printf("       SYM : %10s  OUTFILE => %s %s  INTERFACE %u (min %u)\n", sym.c_str(), pnameTo, addRemoveText, interfaceContacts(), ifaceMin);
            printf("      D_KL :  %-5.3E ( %.3E ) TOTAL ENRGY : %.4E\n", currentKL, lowest_energy, current_energy);
        }

        // update run time parameters
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = totalViolations;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

#ifdef ANNEAL_VERIFY
        {
            std::vector<unsigned int> fullBins(maxbin); unsigned int v;
            calculateModelPrDistributionSym(&bead_indices, &fullBins, workingLimit, v, pModel, pData);
            std::set<unsigned int> tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
            if (fullBins != binCount || pData->getScore(fullBins) != currentKL
                || bruteForceViolations(bead_indices, workingLimit, pModel, violation_limit) != totalViolations
                || calculateNumberOfPairwiseContacts(workingLimit, pModel, tree, bead_indices) != contactsEnergy
                || eulerTour.validateList("refineSym") || eulerTour.getNumberOfComponents() != currentNumberOfComponents
                || std::fabs(rebuildEnergy() - current_energy) > 1e-9*std::fabs(current_energy)){
                std::cout << "VERIFY FAILED (refineSymModel) at step " << numberOfCoolingTempSteps << std::endl; exit(1);
            }
        }
#endif

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        counter++;
    } // end of steps

    std::set<unsigned int> finalSet(bead_indices.begin(), bead_indices.begin() + workingLimit);
    calculateModelParametersSymmetry(&finalSet, pModel);

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
        tempAverageContacts += contacts.contacts(bead_indices[i]);
    }

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
            if (contacts.contacts(bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    logger("PRINTING", "DISTRIBUTION OF CONTACTS");
    for (unsigned int c=1; c<13; c++){
        unsigned int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            if (contacts.contacts(bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }

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

    std::string status = "RAMPING TO TEMP";

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
 * Starting from a filled space, sample down using Add or Remove
 * Starting Model has a connected Euler tour
 *
 * Fixes: status printed with %s from a std::string (undefined behaviour); percentAddRemove (a member) was scaled
 * during the run although this loop never uses it; the starting energy and the energy rebuilt in the ramp omitted
 * beta*violations, which every trial energy includes.
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineSymModelRefine(PointSetModel *pModel, PofRData *pData, std::string nameTo){
    const std::string sym = pModel->getSymmetry();
    const char * pnameTo = nameTo.c_str() ;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.0000001 : (double)highTempStartForCooling;
    srand(time(0));
    std::mt19937 gen(seedForRun(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;
    std::string status = "RAMPING TO TEMP";

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    /*
     * copy model into bead_indices
     */
    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit(); // restrict universe to size of the input model that is going to be refined

    // randomize upto workingLimit and use original
    std::shuffle(bead_indices.begin(), bead_indices.begin()+workingLimit, gen);
    workingLimit = pModel->getBaseWorkingLimit();
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin()+workingLimit);

    /*
     * create set of beads to select from based on all those that are in contact with selected set
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    std::vector<unsigned int> hullPoints(hull.begin(), hull.end());

    if (maxbin < 2){
        throw std::invalid_argument(" INCORRECT MAXBIN NOT INITIALIZED PROPERLY");
    }

    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
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
    auto readyToAddLimit = pModel->getBaseWorkingLimit();

    bool isUpdated = false;
    float acceptRate = 0.5;
    const float inv500slash499 = 499.0f/500.0f, inv500 = 1.0f/500.0f;
    unsigned int original, counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    char addRemoveText[50];

    // incremental bookkeeping
    ContactTracker contacts(pModel->getDirectPointerToNeighborhood(), pModel->getSizeOfNeighborhood(), pModel->getNeighborLimit(), totalBeadsInSphere);
    contacts.reset(bead_indices.begin(), workingLimit);
    SwapSelection selection(bead_indices, workingLimit, totalBeadsInSphere);
    ensureSymmetryCache(pModel);
    ensureClashTable(pModel);
    trkReset(pModel, bead_indices.begin(), workingLimit);
    ifaceArmed = true; ifaceHard = ifaceOn && interfaces.count() >= ifaceMin;   // refinement: hard once the requirement is met
    totalViolations = clashes.count();

    // contacts potential
    double tempContactEnergy;
    double contactsEnergy = contacts.energy();

    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);
    const auto tenPercent = (unsigned int)(0.10*step_limit);

    double contactDistributionWeight = currentKL/contactsEnergy;

    double current_energy = currentKL +
            lambda*connectivityPotential(currentNumberOfComponents) +
            contactDistributionWeight*contactsEnergy +
            beta*totalViolations + ifacePenalty();          // (was missing: trial energies include beta*violations)

    double lowest_energy = current_energy, this_energy;

    logger("alpha", formatNumber(alpha, 3));
    logger("lambda", formatNumber(lambda, 3));
    logger("beta", formatNumber(beta, 3));

    logger("NUMBER OF COMP", std::to_string(currentNumberOfComponents));
    logger("VIOLATIONS", std::to_string(totalViolations));
    logger("D_KL", formatNumber(currentKL, 7));
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

    unsigned int numberOfCoolingTempSteps=0;
    bool doIt=false;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        startTime = std::clock();

        // additional points to expand deadlimit will occur via enlarging CVX Hull
        if (readyToAdd && distribution(gen) < 0.499){ // ADD BEAD?
            original  = hullPoints[randomHull(gen)];

            selection.add(original);                     // swaps to working position, workingLimit += 1
            contacts.add(original);
            addToPrSymAt(workingLimit-1, bead_indices, workingLimit, binCount, pData);
            testKL = pData->getScore(binCount);

            tempContactEnergy = contacts.energy();
            tempNumberOfComponents = eulerTour.addNode(original, pModel);

            trkAdd(original);
            temp_violations = clashes.count();

            this_energy = testKL +
                lambda*connectivityPotential(tempNumberOfComponents) +
                    beta*temp_violations + ifacePenalty() + contactDistributionWeight*tempContactEnergy;

            if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     ADD => %i", 1);
                hull.erase(original);
            } else { // undo changes (rejecting)
                eulerTour.removeNode(original);
                contacts.remove(original);
                trkRemove(original);
                std::sprintf(addRemoveText, "     ADD => %i", 0);
                selection.undoAdd();
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
            }

        } else { // REMOVE BEADS?

            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
            original = bead_indices[selections[0]];
            // test for deletion
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

            removeFromPrSymAt(selection.position(original), bead_indices, workingLimit, binCount, pData);
            selection.remove(original);                  // workingLimit -= 1
            contacts.remove(original); // remove bead for next calculation
            trkRemove(original);

            temp_violations = clashes.count();
            testKL = pData->getScore(binCount);

            tempContactEnergy = contacts.energy();

            this_energy = testKL +
                lambda*connectivityPotential(tempNumberOfComponents) +
                    beta*temp_violations + ifacePenalty() +
                        contactDistributionWeight*tempContactEnergy;

            if (this_energy < current_energy || (std::exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) )) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     REM => %i", original);
                hull.insert(original);
            } else { // undo changes and move to next bead (rejecting)
                contacts.add(original);
                trkAdd(original);
                selection.undoRemove();                  // workingLimit += 1
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                eulerTour.addNode(original, pModel);
            }
        }

        if (numberOfCoolingTempSteps % 1001 == 0){
            printf("      TEMP : %-.3E MAXSTEPS => %.0f (%i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, status.c_str());
            printf("    ACCEPT : %10.5f FAILURES => %6i  TIME : %.5f\n", acceptRate, failures, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
            if (ifaceOn && !ifaceHard && interfaces.count() >= ifaceMin) ifaceHard = true;
            printf("       SYM : %10s  OUTFILE => %s %s  INTERFACE %u (min %u)\n", sym.c_str(), pnameTo, addRemoveText, interfaceContacts(), ifaceMin);
            printf("    E_VIOL :   %.2E VIOL => %d LIMIT : %5i COMP %4i HULL %i  INTERFACE %u (min %u)\n", (beta*totalViolations), totalViolations, workingLimit, currentNumberOfComponents, hullsize, interfaceContacts(), ifaceMin);
            printf("  E_CNTCTS :   %.2E (%.2E)  \n", (contactsEnergy*contactDistributionWeight), contactsEnergy);
            printf("      D_KL :  %-5.3E ( %.3E ) TOTAL ENRGY : %.4E\n", currentKL, lowest_energy, current_energy);
        }

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
            contactsEnergy = tempContactEnergy;
            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            /*
             * do not update hull until target acceptance rate is reached
             */
            if(doIt){
                std::set<unsigned int> inUse(bead_indices.begin(), bead_indices.begin() + workingLimit);
                populateLayeredDeadlimitUsingSet(inUse, hull, pModel);
            }
            hullsize = hull.size();
            hullPoints.assign(hull.begin(), hull.end());
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

            if (workingLimit <= readyToAddLimit){
                readyToAdd=true;
            }
        } else {
            if (!ifaceBlocked){ acceptRate = inv500slash499*acceptRate; }   // a move blocked by the interface constraint is not a failed trial
            if (!ifaceBlocked) failures++;
            ifaceBlocked = false;
        }

        if (!doIt && numberOfCoolingTempSteps > fifteenPercent){
            doIt = true;
        }

        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent ){
            status = "RAMPING";
            contactDistributionWeight = currentKL/contactsEnergy;//*0.733;
            current_energy = currentKL + contactDistributionWeight*contactsEnergy + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + ifacePenalty();
        } else if (numberOfCoolingTempSteps < sixtyfivePercent) {
            status = "EQUILIBRATING";
        }

        if (numberOfCoolingTempSteps > sixtyfivePercent && numberOfCoolingTempSteps%tenPercent == 0) { // reset the weight every 10%
            contactDistributionWeight = currentKL/contactsEnergy*0.5372;
            current_energy = currentKL +
                lambda*connectivityPotential(currentNumberOfComponents) +
                    beta*totalViolations + ifacePenalty() +
                        contactDistributionWeight*contactsEnergy;
        }


        if (numberOfCoolingTempSteps == sixtyfivePercent){
            contactDistributionWeight = currentKL/contactsEnergy*0.372;
            current_energy = currentKL +
                            lambda*connectivityPotential(currentNumberOfComponents) +
                                beta*totalViolations + ifacePenalty() +
                                    contactDistributionWeight*contactsEnergy;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        counter++;
    } // end of steps

    std::set<unsigned int> finalSet(bead_indices.begin(), bead_indices.begin() + workingLimit);
    calculateModelParametersSymmetry(&finalSet, pModel);

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
        tempAverageContacts += contacts.contacts(bead_indices[i]);
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
            if (contacts.contacts(bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (unsigned int c=1; c<13; c++){
        unsigned int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            if (contacts.contacts(bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }

    pData->printKLDivergence(binCount);

    return nameOfModel;
}
