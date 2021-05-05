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
#include "../src/EulerTour/Node.cpp"

class NodeTest : public ::testing::Test {


public:
    Node node1;
    Node node2;
    Node node3;

    NodeTest() : ::testing::Test(),
                    node1(1),
                    node2(2),
                    node3(3){
    }
};


TEST_F(NodeTest, testGetKey){
    ASSERT_EQ(node1.getKey(), 1);
    ASSERT_EQ(node2.getKey(), 2);
}


TEST_F(NodeTest, testAddNeighbor){
    node1.addNeighbor(&node2);
    node1.addNeighbor(&node3); // 1-2-1-3 (1 has two neighbors)

    ASSERT_EQ(node1.getTotalNeighbors(), 2);
    ASSERT_EQ(node2.getTotalNeighbors(), 1);
    ASSERT_EQ(node3.getTotalNeighbors(), 1);
}

TEST_F(NodeTest, testRemoveNeighbor){
    node1.addNeighbor(&node2);
    node1.addNeighbor(&node3); // 1-2-1-3 (1 has two neighbors)

    ASSERT_EQ(node1.getTotalNeighbors(), 2);
    ASSERT_EQ(node2.getTotalNeighbors(), 1);
    ASSERT_EQ(node3.getTotalNeighbors(), 1);

    node3.removeNeighbor(&node1);
    ASSERT_EQ(node3.getTotalNeighbors(), 0);
    ASSERT_EQ(node2.getTotalNeighbors(), 1);
}


TEST_F(NodeTest, isNeighborPresent){
    node1.addNeighbor(&node2);
    node1.addNeighbor(&node3); // 1-2-1-3 (1 has two neighbors)

    ASSERT_TRUE(node1.isNeighborPresent(node2.getKey()));
    ASSERT_TRUE(node1.isNeighborPresent(node3.getKey()));

    node1.removeNeighbor(&node3);
    ASSERT_FALSE(node1.isNeighborPresent(node3.getKey()));
}