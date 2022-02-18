//
// Created by xos81802 on 17/06/2021.
//
// Copyright 2021 Robert P. Rambo
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
#include <sastools/IofQData.h>
#include "gtest/gtest.h"
#include "../src/DensityMapper.cpp"
#include "../src/SHEMapper/Shell.cpp"
#include "../src/SHEMapper/Neighbors.cpp"


class DensityMapperTest : public ::testing::Test {

public:
    std::string filename, iofqdata_filename, filtered, bsa_pdb_model, bsa_test_model;
    // setup
    DensityMapperTest() : ::testing::Test(), filename(  tests::fixture_files("centered_bsa_bead.pdb"), false),
                          iofqdata_filename(tests::fixture_files("BSA_b_refined_sx.dat"), false),
                          filtered(  tests::fixture_files("filtered.pdb"), false),
                          bsa_pdb_model(tests::fixture_files("bsa.pdb"), false),
                          bsa_test_model(tests::fixture_files("centered_SHE_test_model.pdb")){

    }
};


TEST_F(DensityMapperTest, basicConstructorTest){


    float qmax = 0.32;
    float sampling = 1.2;

    DensityMapper dm(filename, qmax, sampling);

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax

    ASSERT_EQ(9, dm.getTotalShells());
}


TEST_F(DensityMapperTest, loadIofQFileConstructorTest){

    float sampling = 2.5;

    IofQData iofqdata_bsa = IofQData(iofqdata_filename, false);
    iofqdata_bsa.extractData();
    iofqdata_bsa.makeWorkingSet();
    auto workingset = iofqdata_bsa.getWorkingSet();

//    for(auto datum : workingset){
//        std::cout << datum.getQ() << " " << datum.getI() << std::endl;
//    }

    float qmax = iofqdata_bsa.getQmax();

    ASSERT_GE(qmax, 0.22);

    DensityMapper dm(filtered, qmax, sampling);

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax
    dm.setBessels(iofqdata_bsa.getQvalues());
    dm.calculateDensityCoefficientAtLMR();

    auto qvalues = iofqdata_bsa.getQvalues();
//    int total_in = qvalues.size();
//    for(int i=0; i<total_in; i++){
//        std::cout << i << " " << qvalues[i] << " " << dm.calculateIntensityAtQ(i) << std::endl;
//    }

   // dm.refineModel(3, 0.05, 10000, workingset);

    ASSERT_EQ(16, dm.getTotalShells());
}


TEST_F(DensityMapperTest, printLatticePointsTest){

    float sampling = 2.5;

    float qmax = 0.22;

    DensityMapper dm(filtered, qmax, sampling);

    dm.printLatticePointsInfo();

}


TEST_F(DensityMapperTest, densityGridTest){

    float sampling = 3.5;

    IofQData iofqdata_bsa = IofQData(iofqdata_filename, false);
    iofqdata_bsa.extractData();
    iofqdata_bsa.makeWorkingSet();
    auto workingset = iofqdata_bsa.getWorkingSet();

    float qmax = iofqdata_bsa.getQmax();

//    DensityMapper dm(filename, qmax, sampling);
    DensityMapper dm(filtered, qmax, sampling);

    // calculate map
    std::uniform_real_distribution<> distribution(0.0,1.0);
    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<float> amps;
    for(int i=0; i<dm.getTotalCenteredCoordinates(); i++){
        amps.push_back((float)distribution(gen));
    }
    dm.setAmplitudes(amps);

    // make map
    dm.createXPLORMap("test");

    ASSERT_EQ(22, dm.getTotalShells());
}



TEST_F(DensityMapperTest, testerTest){

    float sampling = 3.5;

    IofQData iofqdata_bsa = IofQData(iofqdata_filename, false);
    iofqdata_bsa.extractData();
    iofqdata_bsa.makeWorkingSet();
    auto workingset = iofqdata_bsa.getWorkingSet();

    float qmax = iofqdata_bsa.getQmax();

    DensityMapper dm(bsa_test_model, qmax, sampling);

    dm.setBessels(iofqdata_bsa.getQvalues());
    dm.tester(bsa_pdb_model, workingset); //

    // calculate map
    ASSERT_EQ(23, dm.getTotalShells());
}


TEST_F(DensityMapperTest, updateDensityTest){
    float sampling = 2.5;

    IofQData iofqdata_bsa = IofQData(iofqdata_filename, false);
    iofqdata_bsa.extractData();
    iofqdata_bsa.makeWorkingSet();
    auto workingset = iofqdata_bsa.getWorkingSet();

    float qmax = iofqdata_bsa.getQmax();

    DensityMapper dm(bsa_test_model, qmax, sampling);

    std::vector<Datum> & workingSet = const_cast<std::vector<Datum> &>(iofqdata_bsa.getWorkingSet());
    std::vector<Datum> & workingSetSmoothed = const_cast<std::vector<Datum> &>(iofqdata_bsa.getWorkingSetSmoothed());


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(0.0,1.0);

    const unsigned int total_data = workingSet.size();
    Datum * const pWorkingSet = workingSet.data();
    Datum * const pWorkingSetSmoothed = workingSetSmoothed.data();
    std::vector<float> i_calc(total_data);
    std::vector<float> residuals(total_data);
    std::vector<float> qvalues;
    std::vector<float> sigmas_squared;

    for(auto & datum : workingSet){
        qvalues.push_back(datum.getQ());
        sigmas_squared.push_back(1.0f/(datum.getSigma()*datum.getSigma()));
    }

    auto & lattice_points = dm.getLatticePoints();
    const LatticePoint * pLattice = lattice_points.data();
    auto total_lattice_points = lattice_points.size();

    std::vector<float> amplitudes;
    amplitudes.resize(total_lattice_points); // holds the amplitudes in use


    unsigned int total_amplitudes_per_lattice_point = lattice_points[0].getTotalAmplitudes();
    std::vector<unsigned int> selectedIndices(total_lattice_points);
    std::vector<unsigned int> lattice_indices(total_lattice_points);
    std::uniform_int_distribution<unsigned int> randomAmpIndex(0,total_amplitudes_per_lattice_point-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomLatticePt(0,total_lattice_points-1); // guaranteed unbiased

    std::vector<unsigned int> amplitude_indices(total_amplitudes_per_lattice_point);
    for(unsigned int i=0; i<total_amplitudes_per_lattice_point; i++){
        amplitude_indices[i] = i;
    }

    unsigned int amp_index, amp_index2, selected_lattice_pt, selected_lattice_pt2;
    std::vector<float> prior_model_amplitudes(total_lattice_points);

    // set up initial random density model
    float * const pPriors = prior_model_amplitudes.data();
    float * const pAmp = amplitudes.data();

    for(int i=0; i < total_lattice_points; i++){
        pPriors[i] = pLattice[i].getAmplitudeByIndex(randomAmpIndex(gen));
        pAmp[i] = pPriors[i];
        lattice_indices[i] = i;
    }

    dm.createHCPGrid();
    unsigned int long totalKept = dm.getTotalKeptHCPSamplingPoints();
    std::vector<float> hcp_electron_densities(totalKept);


    // use aligned memory for squared amplitudes and debye_factors
    auto * hcp_squared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalKept*(totalKept-1)/2 + totalKept)), 16);
    auto * testsquared_amplitudes = (float *)_aligned_malloc(sizeof(float)*((totalKept*(totalKept-1)/2 + totalKept)), 16);

    //takes less than 100 usec
    float ave = dm.populateDensities(prior_model_amplitudes, hcp_electron_densities);
    dm.populateSquaredAmplitudes(totalKept, hcp_electron_densities, hcp_squared_amplitudes);
    float ave2 = dm.populateDensities(prior_model_amplitudes, hcp_electron_densities);

    std::clock_t startTime;
    double runtime;

    // make a change in prior and check that
    /*
     *
     */
    std::vector<unsigned int> indicesToModify(1);
    for(int i=0; i<1000; i++){
        selected_lattice_pt = randomLatticePt(gen);
        indicesToModify[0] = selected_lattice_pt;

        std::shuffle(amplitude_indices.begin(), amplitude_indices.end(), gen);
        amp_index = amplitude_indices[0];
//
        float * pPrior = &pPriors[selected_lattice_pt];
        *pPrior = pLattice[selected_lattice_pt].getAmplitudeByIndex(amp_index); // set new amplitudes

        if ( fabs(pAmp[selected_lattice_pt] - *pPrior) < 0.0001f){ // make sure we pick a different value
            amp_index = amplitude_indices[1];
            *pPrior = pLattice[selected_lattice_pt].getAmplitudeByIndex(amp_index);
        }

        ave = dm.populateDensities(prior_model_amplitudes, hcp_electron_densities);
        dm.populateSquaredAmplitudes(totalKept, hcp_electron_densities, testsquared_amplitudes);
        dm.updateDensities(indicesToModify, prior_model_amplitudes, hcp_electron_densities, testsquared_amplitudes);
//
        for(unsigned int t=0; t < (totalKept*(totalKept-1)/2 + totalKept); t++){
            std::string message = "Round " + std::to_string(i) + " " + std::to_string(t);
            ASSERT_FLOAT_EQ(hcp_squared_amplitudes[t], testsquared_amplitudes[t]) << message;
        }

        // undo
        *pPrior = pAmp[selected_lattice_pt]; // reverse change in prior_model_amplitude
        dm.populateSquaredAmplitudes(totalKept, hcp_electron_densities, testsquared_amplitudes);
       dm.updateDensities(indicesToModify, prior_model_amplitudes, hcp_electron_densities, testsquared_amplitudes);
    }

    std::cout << "finished " << std::endl;
    _aligned_free(hcp_squared_amplitudes);
    _aligned_free(testsquared_amplitudes);

}

TEST_F(DensityMapperTest, checkNeighbors){
    float sampling = 2.5;
    float qmax = 0.22;

    DensityMapper dm(filtered, qmax, sampling);

    // lattice_points should be filled along with neighbors

    std::set<unsigned int> indices;

    std::vector<LatticePoint> & pLats = dm.getLatticePoints();
    std::vector<LatticePoint> omp_lattice_points(dm.getLatticePoints());



    unsigned int total = pLats.size();
    for(unsigned int i=0; i<total; i++){
        indices.insert(i);
    }

    ASSERT_FALSE(dm.checkConnectivity(indices, omp_lattice_points));

    // make connected set
    std::set<unsigned int> connected;

    for(auto & lat : pLats){
        if (lat.getTotalNeighbors() > 0 && lat.getIndex() != 273 && lat.getIndex() != 277){
            connected.insert(lat.getIndex());
        }
    }

    ASSERT_TRUE(dm.checkConnectivity(connected, omp_lattice_points));

}

//TEST_F(DensityMapperTest, createHCPGrid){
//
//    float sampling = 2.35;
//
//    IofQData iofqdata_bsa = IofQData(iofqdata_filename, false);
//    iofqdata_bsa.extractData();
//    iofqdata_bsa.makeWorkingSet();
//    auto workingset = iofqdata_bsa.getWorkingSet();
//
//    float qmax = iofqdata_bsa.getQmax();
//
//    DensityMapper dm(bsa_test_model, qmax, sampling);
//
//    dm.setBessels(iofqdata_bsa.getQvalues());
//    dm.createHCPGrid(); //
//    dm.setYlms();
//
//    // calculate map
//    ASSERT_EQ(23, dm.getTotalShells());
//}