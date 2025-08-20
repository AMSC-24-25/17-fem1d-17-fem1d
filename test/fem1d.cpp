#include <gtest/gtest.h>
#include "fem1d.hpp"
#include <cmath>
#include <vector>

class FEM1DTest : public ::testing::Test {
protected:
    // Helper function to compute L2 error
    double computeL2Error(const std::vector<double>& numerical, 
                          const std::vector<double>& mesh,
                          std::function<double(double)> exact) {
        double error = 0.0;
        for (size_t i = 0; i < numerical.size(); i++) {
            double diff = numerical[i] - exact(mesh[i]);
            error += diff * diff;
        }
        return std::sqrt(error / numerical.size());
    }
};

TEST_F(FEM1DTest, Test_test) {
    GTEST_SUCCEED(); 
}

