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
#include "gtest/gtest.h"
#include "../src/SubUnit.cpp"

TEST(SubUnitTest, testGetKey){
    std::random_device rd;
    std::mt19937 gen(rd());

    float radius = 3.75f;
    PointSetModel basic(150, radius);

    unsigned int subUnitWorkingLimit = 100;
    unsigned int totalBeadsInSphere = basic.getTotalNumberOfBeadsInUniverse();

    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomInUniverse (0, totalBeadsInSphere-1);

    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    unsigned int * ptr = (totalBeadsInSphere != 0) ? &subUnit_indices.front() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit);
    unsigned int translateTo = subUnit_indices[randomInUniverse(gen)];
    SubUnit sub1(0, subUnit_indices[randomIndex(gen)],translateTo, subUnitWorkingLimit,subUnit_indices, &basic);

    ASSERT_EQ(sub1.getSubUnitIndex(), 0);
}

TEST(SubUnitTest, testCopyAssignment){
    std::random_device rd;
    std::mt19937 gen(rd());

    float radius = 3.75f;
    PointSetModel basic(150, radius);

    unsigned int subUnitWorkingLimit = 100;
    unsigned int totalBeadsInSphere = basic.getTotalNumberOfBeadsInUniverse();

    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomInUniverse (0, totalBeadsInSphere-1);

    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    unsigned int * ptr = (totalBeadsInSphere != 0) ? &subUnit_indices.front() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin()+subUnitWorkingLimit);
    unsigned int translateTo = subUnit_indices[randomInUniverse(gen)];
    SubUnit sub1(0, subUnit_indices[randomIndex(gen)],translateTo, subUnitWorkingLimit,subUnit_indices, &basic);

    SubUnit sub2(sub1);

    ASSERT_EQ(sub1.getSubUnitIndex(), sub2.getSubUnitIndex());
    ASSERT_EQ(sub2.getTranslateToIndex(), sub1.getTranslateToIndex());

    SubUnit sub3 = sub1;

    ASSERT_EQ(sub1.getSubUnitIndex(), sub3.getSubUnitIndex());
    ASSERT_EQ(sub3.getTranslateToIndex(), sub1.getTranslateToIndex());

}