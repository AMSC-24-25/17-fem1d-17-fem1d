#include <gtest/gtest.h>
#include "function.hpp"
#include "vector.hpp"

class VectorTest : public ::testing::Test {};

TEST_F(VectorTest, ElementwiseOpsAndSize) {
	using F = Function<1>;
	// v(x) = [x, x^2]
	Vector<1> v({
		F([](Point<1> x){ return x[0]; }, [](Point<1>){ return 0.0; }),
		F([](Point<1> x){ return x[0]*x[0]; }, [](Point<1> x){ return 2.0*x[0]; })
	});
	// w(x) = [1, 2x]
	Vector<1> w({
		F([](Point<1>){ return 1.0; }, [](Point<1>){ return 0.0; }),
		F([](Point<1> x){ return 2.0*x[0]; }, [](Point<1>){ return 2.0; })
	});

	EXPECT_EQ(v.size(), 2u);
	EXPECT_EQ(w.size(), 2u);

	Point<1> p(3.0);

	// Addition
	auto a = v + w; // [x+1, x^2+2x]
	EXPECT_DOUBLE_EQ(a[0](p), 3.0 + 1.0);
	EXPECT_DOUBLE_EQ(a[1](p), 9.0 + 6.0);
	EXPECT_DOUBLE_EQ(a[0].dx_value(p), 0.0);
	EXPECT_DOUBLE_EQ(a[1].dx_value(p), 2.0 * 3.0 + 2.0);

	// Subtraction (implemented via + and scalar *)
	auto s = v - w; // [x-1, x^2-2x]
	EXPECT_DOUBLE_EQ(s[0](p), 3.0 - 1.0);
	EXPECT_DOUBLE_EQ(s[1](p), 9.0 - 6.0);

	// Elementwise product
	auto m = v * w; // [x*1, x^2*2x] = [x, 2x^3]
	EXPECT_DOUBLE_EQ(m(p), 3.0 + 2.0 * 27.0);
}

TEST_F(VectorTest, OutOfRangeAccessThrows) {
	using F = Function<1>;
	Vector<1> v({ F([](Point<1>){ return 0.0; }, [](Point<1>){ return 0.0; }) });
	EXPECT_THROW((void)v[1], std::out_of_range);
}

TEST_F(VectorTest, SizeMismatchOpsThrow) {
	using F = Function<1>;
	Vector<1> v({ F([](Point<1>){ return 0.0; }, [](Point<1>){ return 0.0; }) });
	Vector<1> w({
		F([](Point<1>){ return 1.0; }, [](Point<1>){ return 0.0; }),
		F([](Point<1> x){ return x[0]; }, [](Point<1>){ return 0.0; })
	});

	EXPECT_THROW((void)(v + w), std::invalid_argument);
	EXPECT_THROW((void)(v - w), std::invalid_argument);
}

