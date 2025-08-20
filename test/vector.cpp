#include <gtest/gtest.h>
#include <cmath>
#include <vector>

TEST(VectorTest, Dotproduct) {
    std::vector<double> v1{1.0, 2.0, 3.0};
    std::vector<double> v2{4.0, 5.0, 6.0};
    
    double dot = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        dot += v1[i] * v2[i];
    }
    
    EXPECT_NEAR(dot, 32.0, 1e-12);
}

TEST(VectorTest, Norm) {
    std::vector<double> v{3.0, 4.0};
    
    double norm = 0.0;
    for (const auto& val : v) {
        norm += val * val;
    }
    norm = sqrt(norm);
    
    EXPECT_NEAR(norm, 5.0, 1e-12);
}

