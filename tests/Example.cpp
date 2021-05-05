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
#include "../src/utils/Example.cpp"

TEST(WriteExampleTest,checkWriteAnchorFile){
    Example example("anchor");
    // need working directory
    // file is written to CMAKE_RUNTIME_OUTPUT_DIRECTORY
    std::string output = tests::fixture("anchor_example.txt");
    FILE* f = fopen(output.c_str(), "r");
    ASSERT_TRUE(f != nullptr);
}

TEST(WriteExampleTest,checkWriteHelicalFile){
    Example example("helical");
    // need working directory
    // file is written to CMAKE_RUNTIME_OUTPUT_DIRECTORY
    std::string output = tests::fixture("helical_params.txt");
    FILE* f = fopen(output.c_str(), "r");
    ASSERT_TRUE(f != nullptr);
}