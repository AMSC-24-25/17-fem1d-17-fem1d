//
// Created by federico on 06/12/24.
//

#include "gtest/gtest.h"
#include "boundary_cond.hpp"
#include "function.hpp"

class BoundaryCondTest : public ::testing::Test {
protected:
    // Set up the Calculator object before each test
    BoundaryCond cond;

    BoundaryCondTest() : cond(true, ZeroFunction()) {}

};

TEST_F(BoundaryCondTest, ConstructorTest) {
    EXPECT_EQ(cond.isNeumann(), true);
}
