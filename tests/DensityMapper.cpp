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