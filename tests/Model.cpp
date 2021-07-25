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
#include "../src/PointSetModel.cpp"

class ModelTest : public ::testing::Test {


public:
    PointSetModel standardRectangularModel;

    ModelTest() : ::testing::Test(),
                  standardRectangularModel(4.13, 50, 60, 65){
    }
};



TEST_F(ModelTest, testGetKey){
    ASSERT_EQ(standardRectangularModel.getZaxis(), 65);
    ASSERT_EQ(standardRectangularModel.getYaxis(), 60);
    ASSERT_EQ(standardRectangularModel.getXaxis(), 50);
}

TEST_F(ModelTest, getNeighborCutoff){
    ASSERT_NEAR(standardRectangularModel.getNeighborCutOffLimit(), 2.001f*4.13, 0.01);
}



TEST_F(ModelTest, getBeadVolume){
    float radius = 4.13f;
    ASSERT_NEAR(standardRectangularModel.getBeadVolume(), (4.0f/3.0f*radius*radius*radius*M_PI), 0.01); // 4/3*PI*r^3)
}

TEST_F(ModelTest, getSymmetry){
    std::string text = "C1";
    ASSERT_STREQ(standardRectangularModel.getSymmetry().c_str(), text.c_str());
}

TEST_F(ModelTest, checkSymmetrySetter){
    std::string cstring = "C";

    standardRectangularModel.setSymmetryParameters("C2");
    ASSERT_STREQ(standardRectangularModel.getSymmetryGroup().c_str(), cstring.c_str());
    ASSERT_EQ(standardRectangularModel.getSymmetryIndex(), 2);
    ASSERT_EQ(standardRectangularModel.getNumberOfSubUnits(), 2);

    auto totalSubUnits = standardRectangularModel.getNumberOfSubUnits();
    auto pCos = standardRectangularModel.getSubUnitAnglesCos();
    auto pSin = standardRectangularModel.getSubUnitAnglesSin();

    // check angles
    double angle = M_PI*2.0d/(double)standardRectangularModel.getSymmetryIndex();
    ASSERT_EQ(pCos[0], 1.0f);
    ASSERT_EQ(pSin[0], 0.0f);
    ASSERT_EQ(pCos[1], (float)std::cos(angle));
    ASSERT_EQ(pSin[1], (float)std::sin(angle));


    PointSetModel d2Model(4.13, 150, 150, 165);
    cstring = "D";
    d2Model.setSymmetryParameters("D2");
    ASSERT_STREQ(d2Model.getSymmetryGroup().c_str(), cstring.c_str());
    ASSERT_EQ(d2Model.getSymmetryIndex(), 2);
    ASSERT_EQ(d2Model.getNumberOfSubUnits(), 4);
}

TEST_F(ModelTest, checkBasicConstructor){
    float radius = 3.75f;
    PointSetModel basic(150, radius);
    std::string text = "C1";
    ASSERT_STREQ(basic.getSymmetry().c_str(), text.c_str());
    ASSERT_NEAR(basic.getBeadVolume(), (4.0f/3.0f*radius*radius*radius*M_PI), 0.01); // 4/3*PI*r^3)

    unsigned int number_of_beads = basic.getTotalNumberOfBeadsInUniverse();
    unsigned long int totalDistances = ( ((unsigned long int)number_of_beads*(number_of_beads-1)) / 2);
    ASSERT_EQ(basic.getTotalDistances(), totalDistances);
}

TEST_F(ModelTest, checkConvertToBin){

    PofRData ribo( tests::fixture_files("30S_0p22_pr.dat"), false);
    float bead_radius =  (float)(ribo.getBinWidth()/2.0f);
    PointSetModel model(ribo.getDmax()*1.2, bead_radius);
    // pModel->populateBins(pData)
    unsigned short int maxbin = model.populateBins(&ribo);
    ASSERT_EQ(maxbin, 22); // largest bin based on data and size of universe;
}