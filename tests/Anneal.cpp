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
#include <support.hpp>
#include "gtest/gtest.h"
#include "../src/Anneal.cpp"
#include "../src/Component.cpp"
#include "../src/EulerTour/EulerTour.cpp"

class AnnealerTest : public ::testing::Test {

public:
    PofRData ribo30S;
    // setup
    AnnealerTest() : ::testing::Test(), ribo30S( PofRData( tests::fixture_files("30S_0p22_pr.dat"), false)) {

    }
};

TEST_F(AnnealerTest, basicPrLoadTest){
    float radius = (float)(ribo30S.getBinWidth()/2.0);
    float interconnectivityCutOff = radius * 2.001f;
    PointSetModel basic(ribo30S.getDmax()*1.2f, radius);
    std::string text ="tmp";
    Anneal anneal(0.001, 0.44, 1234, text, 0.01, 0.01, 0.02, 0.03, 0.02, 45, 0.12, interconnectivityCutOff);
    anneal.fillPrBinsAndAssignTotalBin(&basic, &ribo30S);

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax
    ASSERT_EQ(anneal.getMaxBin(), 23);
}

TEST_F(AnnealerTest, getUseableNeighborFromSetTest){

    float radius = (float)(ribo30S.getBinWidth()/2.0);
    float interconnectivityCutOff = radius * 2.001f;
    PointSetModel basic(ribo30S.getDmax()*1.2f, radius);

    const unsigned int totalBeadsInSphere = basic.getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    unsigned int workingLimit = (unsigned int)std::floor(totalBeadsInSphere*0.2);

    /*
     * create initial model
     * randomize bead indices
     * sort to workingLimit
     */
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads;
    for(unsigned int i=0; i<workingLimit; i++){
        beads.insert(ptr[i]);
    }

    std::string text ="tmp";
    Anneal anneal(0.001, 0.44, 1234, text, 0.01, 0.01, 0.02, 0.03, 0.02, 45, 0.12, interconnectivityCutOff);

    for(unsigned int i=0; i<100; i++){
        unsigned int index = anneal.getUseableNeighborFromSet(&beads, &basic, bead_indices[i]);
        ASSERT_NE(index, basic.getNeighborLimit());
    }

}


TEST_F(AnnealerTest, connectivityPotentialTest){

    float radius = (float)(ribo30S.getBinWidth()/2.0);
    float interconnectivityCutOff = radius * 2.001f;
    PointSetModel basic(ribo30S.getDmax()*1.2f, radius);

    std::string text ="tmp";
    Anneal anneal(0.001, 0.44, 1234, text, 0.01, 0.01, 0.02, 0.03, 0.02, 45, 0.12, interconnectivityCutOff);

    ASSERT_LE(anneal.connectivityPotential(1), anneal.connectivityPotential(2));
    ASSERT_LE(anneal.connectivityPotential(2), anneal.connectivityPotential(3));
    ASSERT_LE(anneal.connectivityPotential(3), anneal.connectivityPotential(4));
    ASSERT_LE(anneal.connectivityPotential(4), anneal.connectivityPotential(5));

}



TEST_F(AnnealerTest, updateASATempTest){

    float radius = (float)(ribo30S.getBinWidth()/2.0);
    float interconnectivityCutOff = radius * 2.001f;
    PointSetModel basic(ribo30S.getDmax()*1.2f, radius);

    const unsigned int totalBeadsInSphere = basic.getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    unsigned int workingLimit = (unsigned int)std::floor(totalBeadsInSphere*0.2);

    /*
     * create initial model
     * randomize bead indices
     * sort to workingLimit
     */
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads;
    for(unsigned int i=0; i<workingLimit; i++){
        beads.insert(ptr[i]);
    }

    std::string text ="tmp";
    Anneal anneal(0.001, 0.44, 1234, text, 0.01, 0.01, 0.02, 0.03, 0.02, 45, 0.12, interconnectivityCutOff);
    unsigned int index = 16;
    float evalMax = 100.0f;
    double temp=32, inv_temp;

    double oldtemp = temp;

    anneal.updateASATemp(index, evalMax, 0.44, temp, inv_temp);
    ASSERT_LE(temp, oldtemp);

    index = 66;
    oldtemp = temp;
    anneal.updateASATemp(index, evalMax, 0.44, temp, inv_temp);
    ASSERT_LE(temp, oldtemp);
}



TEST_F(AnnealerTest, removeFromContactsDistribution){

    float radius = (float)(ribo30S.getBinWidth()/2.0);
    float interconnectivityCutOff = radius * 2.001f;
    PointSetModel basic(ribo30S.getDmax()*1.2f, radius);

    const unsigned int totalBeadsInSphere = basic.getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    unsigned int workingLimit = (unsigned int)std::floor(totalBeadsInSphere*0.2);

    /*
     * create initial model
     * randomize bead indices
     * sort to workingLimit
     */
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads_in_use_tree;
    for(unsigned int i=0; i<workingLimit; i++){
        beads_in_use_tree.insert(ptr[i]);
    }

    std::string text ="tmp";
    Anneal anneal(0.001, 0.44, 1234, text, 0.01, 0.01, 0.02, 0.03, 0.02, 45, 0.12, interconnectivityCutOff);

    unsigned int original;
    std::vector<double> controlContactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);
    anneal.populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, &basic);

    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    for(unsigned int i=0; i<1000; i++){

        original = bead_indices[randomIndex(gen)];

        anneal.removeFromContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, &basic);

        beads_in_use_tree.erase(original);
        anneal.populateContactsDistribution(controlContactsDistributionOfModel, &beads_in_use_tree, &basic);
        for(unsigned int i=0; i<13; i++){
            ASSERT_EQ(tempContactsDistributionOfModel[i],controlContactsDistributionOfModel[i]) << " REMOVE " << i;
        }
        //undo
        beads_in_use_tree.insert(original);
        anneal.addToContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, &basic);
        anneal.populateContactsDistribution(controlContactsDistributionOfModel, &beads_in_use_tree, &basic);

        for(unsigned int i=0; i<13; i++){
            ASSERT_EQ(tempContactsDistributionOfModel[i],controlContactsDistributionOfModel[i]) << " ADD " << i;
        }
    }

}

TEST_F(AnnealerTest, estimateMagnitudeOfDifferenceContactsPotential){

    float radius = (float)(ribo30S.getBinWidth()/2.0);
    float interconnectivityCutOff = radius * 2.001f;
    PointSetModel basic(ribo30S.getDmax()*1.2f, radius);

    const unsigned int totalBeadsInSphere = basic.getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    unsigned int workingLimit = (unsigned int)std::floor(totalBeadsInSphere*0.2);

    /*
     * create initial model
     * randomize bead indices
     * sort to workingLimit
     */
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads_in_use_tree;
    for(unsigned int i=0; i<workingLimit; i++){
        beads_in_use_tree.insert(ptr[i]);
    }

    std::string text ="tmp";
    Anneal anneal(0.001, 0.44, 1234, text, 0.01, 0.01, 0.02, 0.03, 0.02, 45, 0.12, interconnectivityCutOff);

    unsigned int original;
    std::vector<double> controlContactsDistributionOfModel(13);
    anneal.populateContactsDistribution(controlContactsDistributionOfModel, &beads_in_use_tree, &basic);


    std::vector<unsigned int> binCount(100);        // smallish vector, typically < 50
    anneal.calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, &basic);

    double delta = anneal.estimateMagnitudeOfDifferenceContactsPotential(
            workingLimit, 
            bead_indices,
            binCount,
            &basic);
    
    std::cout << " delta " << delta << std::endl;
    ASSERT_NEAR(delta, 0.0d, 0.00001);
    //ASSERT_GT(delta, 0);
}