#include <gtest/gtest.h>
#include "fem.hpp"
#include "grid1D.hpp"
#include "function.hpp"
#include "boundary_conditions.hpp"
#include "quadrature.hpp"
#include "point.hpp"

#include <cmath>
#include <fstream>

// Helper: exact solution u(x) = x(1-x) on [0,1], so -u'' = 2.
// We test with diffusion = 1, transport = 0, reaction = 0, Dirichlet u(0)=u(1)=0.
static double exact_u(double x) { return x*(1.0 - x); }

TEST(Fem1D, Poisson_Dirichlet_Uxx_2) {
    const int N = 41;              // number of nodes
    const double a = 0.0, b = 1.0; // domain [a,b]
    const double h = (b - a) / (N - 1);

    // Build a uniform 1D grid and convert to Grid<1>
    Grid1D builder(a, b, N);
    Grid<1> grid = static_cast<Grid<1>>(builder);

    // PDE coefficients and RHS
    // Forcing f(x) = 2
    Function<1,1> forcing([](Point<1> p) { return 2.0; });
    // Diffusion D(x) = 1
    Function<1,1> diffusion([](Point<1> p) { return 1.0; });
    // Transport b(x) = 0
    Function<1,1> transport([](Point<1> p) { return 0.0; });
    // Reaction r(x) = 0
    Function<1,1> reaction([](Point<1> p) { return 0.0; });

    // Boundary conditions: Dirichlet u=0 on both endpoints (tags 0 and 1 are used by Grid1D)
    BoundaryConditions<1,1> bc;
    bc.addDirichlet(0, 0.0);
    bc.addDirichlet(1, 0.0);

    // Quadrature: order-2 should be enough for linear elements
    QuadratureRule<1> quadrature = OrderTwoQuadrature<1>();

    // Build FEM object
    Fem<1> fem(grid, forcing, diffusion, transport, reaction, bc, quadrature);

    // Assemble and solve
    EXPECT_NO_THROW(fem.assemble());
    EXPECT_NO_THROW(fem.solve());

    // Fetch solution
    const Eigen::VectorXd& u = fem.getSolution();
    ASSERT_EQ(u.size(), N) << "Solution size should equal number of nodes.";

    // Check BCs are enforced (tolerant threshold)
    EXPECT_NEAR(u[0], 0.0, 1e-10);
    EXPECT_NEAR(u[N-1], 0.0, 1e-10);

    // Compute max nodal error against the exact solution
    double max_err = 0.0;
    for (int i = 0; i < N; ++i) {
        double x = a + i * h;
        max_err = std::max(max_err, std::abs(u[i] - exact_u(x)));
    }

    // For linear elements in 1D with a modest mesh, expect ~O(h^2) accuracy.
    // With N=41, h=0.025, O(h^2) ~ 6.25e-4. We allow some slack.
    EXPECT_LT(max_err, 5e-3) << "Max nodal error too large for Poisson 1D with Dirichlet BCs.";
}

TEST(Fem1D, OutputCsv_WritesFile) {
    const int N = 11;
    Grid1D builder(0.0, 1.0, N);
    Grid<1> grid = static_cast<Grid<1>>(builder);

    Function<1,1> forcing([](Point<1> p) { return 2.0; });
    Function<1,1> diffusion([](Point<1> p) { return 1.0; });
    Function<1,1> transport([](Point<1> p) { return 0.0; });
    Function<1,1> reaction([](Point<1> p) { return 0.0; });

    BoundaryConditions<1,1> bc;
    bc.addDirichlet(0, 0.0);
    bc.addDirichlet(1, 0.0);

    QuadratureRule<1> quadrature = OrderTwoQuadrature<1>();

    Fem<1> fem(grid, forcing, diffusion, transport, reaction, bc, quadrature);
    fem.assemble();
    fem.solve();

    // Write CSV in the build/test_output directory (CMake usually runs tests from a build dir)
    const std::string outname = "fem_test_output.csv";
    EXPECT_NO_THROW(fem.outputCsv(outname));

    // Verify the file exists and has content
    std::ifstream fin(outname);
    ASSERT_TRUE(fin.good()) << "CSV file was not created.";
    std::string line;
    int line_count = 0;
    while (std::getline(fin, line)) ++line_count;

    // Expect at least N lines (header + N nodes). Be lenient due to possible header variations.
    EXPECT_GE(line_count, N) << "CSV should contain at least one row per node.";
}
