//
// Created by xos81802 on 27/06/2021.
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
#include <sastools/include/IofQData.h>
#include "gtest/gtest.h"
#include "../src/SHEMapper/LatticePoint.cpp"


class LatticePointTest : public ::testing::Test {

public:
    std::string filename;
    // setup
    LatticePointTest() : ::testing::Test(), filename(  tests::fixture_files("centered_bsa_bead.pdb"), false) {

    }
};


TEST_F(LatticePointTest, basicConstructorTest){

    LatticePoint lp0(0,5);
    LatticePoint lp1(1,5);
    LatticePoint lp2(2,5);

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax

    ASSERT_EQ(0, lp0.getIndex());
    ASSERT_EQ(1, lp1.getIndex());
    ASSERT_EQ(2, lp2.getIndex());
}


TEST_F(LatticePointTest, basicCounterTest){

    LatticePoint lp0(0,5);

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax

    ASSERT_EQ(0, lp0.getIndex());
}


TEST_F(LatticePointTest, basicProbTest){

    LatticePoint lp0(0,5);

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax

    ASSERT_GE(lp0.getProbability(0), 0);
    ASSERT_GE(lp0.getProbability(1), 0);
    ASSERT_GE(lp0.getProbability(2), 0);
}

TEST_F(LatticePointTest, sumProbabilityTest){

    LatticePoint lp0(0,5);

    std::uniform_real_distribution<> distribution(0.0,1.0);
    std::random_device rd;
    std::mt19937 gen(rd());

    // maxBin should be adjusted to 19 for ribo data, 150 is too small for dmax

    float sum=0;
    for(int i=0; i<(6); i++){
        sum += lp0.getProbability(i);
    }

    ASSERT_LE(sum, 1);
}

TEST_F(LatticePointTest, addToCounterTest){

    LatticePoint lp0(0,5);

    std::uniform_real_distribution<> distribution(0.0,1.0);
    std::random_device rd;
    std::mt19937 gen(rd());

    for(int i=0; i<100; i++){
        lp0.addToCounter(distribution(gen));
    }

    lp0.updateProbabilities();

    float sum=0;
    for(int i=0; i<(6); i++){
        sum += lp0.getProbability(i);
    }

    ASSERT_NEAR(sum, 1, 0.01);
}


TEST_F(LatticePointTest, testCDFTest){

    /*
     * i AMP    probability
     * 0 0.0 0.0    < v <= 0.1667
     * 1 0.2 0.1667 < v <= 0.3333
     * 2 0.4 0.3333 < v <= 0.5
     * 3 0.6 0.5    < v <= 0.6666
     * 4 0.8 0.6667 < v <= 0.833
     * 5 1.0 0.8333 < v <= 1.0
     */
    LatticePoint lp0(0,5);
    ASSERT_EQ(lp0.CDF(0.1), 0);
    ASSERT_NEAR(lp0.CDF(0.2), 0.2, 0.01);
    ASSERT_NEAR(lp0.CDF(0.33), 0.2, 0.01);
    ASSERT_NEAR(lp0.CDF(0.49), 0.4, 0.01);
    ASSERT_NEAR(lp0.CDF(0.65), 0.6, 0.01);
    ASSERT_NEAR(lp0.CDF(0.83), 0.8, 0.01);
    ASSERT_NEAR(lp0.CDF(0.85), 1.0, 0.01);
}