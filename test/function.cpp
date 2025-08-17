#include <gtest/gtest.h>
#include "function.hpp"

class functionTest : public ::testing::Test {
};

TEST_F(functionTest, functionTest_T1_Test) {
    Function<1> f1([](Point<1> x) { return x[0]; }, [](Point<1> x) { return 0; });
    Function<1> f2([](Point<1> x) { return x[0] * x[0]; }, [](Point<1> x) { return 2 * x[0]; });

    EXPECT_EQ(f1(2.0), 2.0);
    EXPECT_EQ(f1.dx_value(2.0), 0.0);
    EXPECT_EQ(f2(2.0), 4.0);
    EXPECT_EQ(f2.dx_value(2.0), 4.0);
}


TEST_F(functionTest, functionTest_T2_Test) {
    // Create function objects
    Function<1> f1([](Point<1> x) { return std::sin(x[0]); }, [](Point<1> x) { return std::cos(x[0]); });
    Function<1> f2([](Point<1> x) { return std::exp(x[0]); }, [](Point<1> x) { return std::exp(x[0]); });

    // Test function evaluation
    EXPECT_NEAR(f1(0.0), 0.0, 1e-10);
    EXPECT_NEAR(f1(M_PI/2), 1.0, 1e-10);
    EXPECT_NEAR(f1.dx_value(0.0), 1.0, 1e-10);

    // Test second function
    EXPECT_NEAR(f2(0.0), 1.0, 1e-10);
    EXPECT_NEAR(f2(1.0), std::exp(1.0), 1e-10);
    EXPECT_NEAR(f2.dx_value(0.0), 1.0, 1e-10);
}

// TEST_F(functionTest, Operators_Sum_Product_Functions_And_Scalars) {
//     // f(x) = x, f'(x) = 0 (as set in original test)
//     Function<1> f([](Point<1> x) { return x[0]; }, [](Point<1> /*x*/) { return 0.0; });
//     // g(x) = x^2, g'(x) = 2x
//     Function<1> g([](Point<1> x) { return x[0] * x[0]; }, [](Point<1> x) { return 2.0 * x[0]; });

//     Point<1> p(2.0);

//     // Sum h = f + g => h(x) = x + x^2, h'(x) = f' + g' = 0 + 2x
//     auto h = f + g;
//     EXPECT_DOUBLE_EQ(h(p), 2.0 + 4.0);
//     EXPECT_DOUBLE_EQ(h.dx_value(p), 4.0);

//     // Product q = f * g => q(x) = x^3, q'(x) = f*g' + f'*g = x*(2x) + 0 = 2x^2
//     auto q = f * g;
//     EXPECT_DOUBLE_EQ(q(p), 8.0);
//     EXPECT_DOUBLE_EQ(q.dx_value(p), 2.0 * 4.0);

//     // Scalar sum r = f + c, r'(x) = f'(x)
//     double c = 3.0;
//     auto r = f + c;
//     EXPECT_DOUBLE_EQ(r(p), 2.0 + 3.0);
//     EXPECT_DOUBLE_EQ(r.dx_value(p), 0.0);

//     // Scalar product s = c * g, s'(x) = c * g'(x)
//     auto s = g * c;
//     EXPECT_DOUBLE_EQ(s(p), 3.0 * 4.0);
//     EXPECT_DOUBLE_EQ(s.dx_value(p), 3.0 * 4.0);
// }