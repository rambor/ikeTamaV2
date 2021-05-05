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
#include "../src/Objective.cpp"

TEST(ObjectiveTests, addDatasetTest){
    Objective obj;
    obj.addDataObject(tests::fixture_files("30S_0p22_pr.dat"));
    auto pData = obj.getMainDataset();

    ASSERT_TRUE(pData->getIsPr());
    PofRData * pD = dynamic_cast<PofRData*>(obj.getMainDataset());
    ASSERT_NEAR(pD->getRave(), 85.892, 0.001);
}

TEST(ObjectiveTests, objectiveVectorTest){

}