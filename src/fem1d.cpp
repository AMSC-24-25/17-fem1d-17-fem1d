#include "fem1d.hpp"
#include <Eigen/Dense>
#include <vector>

void gauss_legendre_quadrature(int nQuadraturePoints, std::vector<double>& quad_points, std::vector<double>& quad_weights) {
    quad_points.clear();
    quad_weights.clear();

    // Precompute sqrt values for accuracy
    const double sqrt3 = std::sqrt(3.0);
    const double sqrt3_5 = std::sqrt(3.0/5.0);

    if (nQuadraturePoints == 1) {
        quad_points = {0.5};
        quad_weights = {1.0};
    } else if (nQuadraturePoints == 2) {
        quad_points = {(1.0 - 1.0/sqrt3)/2.0, (1.0 + 1.0/sqrt3)/2.0};
        quad_weights = {0.5, 0.5};
    } else if (nQuadraturePoints == 3) {
        quad_points = {(1.0 - sqrt3_5)/2.0, 0.5, (1.0 + sqrt3_5)/2.0};
        quad_weights = {5.0/18.0, 4.0/9.0, 5.0/18.0};
    } else {
        throw std::runtime_error("nQuadraturePoints must be 1, 2, 3, or 4 for hardcoded quadrature.");
    }
}
void Fem1D::assemble() {
    std::vector<Triplet> triplets;
    triplets.reserve(3 * mesh.getN() - 2);

    unsigned int nQuadraturePoints = 2;
    std::vector<double> quad_points(nQuadraturePoints);
    std::vector<double> quad_weights(nQuadraturePoints);
    gauss_legendre_quadrature(nQuadraturePoints, quad_points, quad_weights);

    for(int elem=0; elem < mesh.getN()-1; ++elem) {
        MatrixXd mat(2, 2); 
        mat.setZero();

        VectorXd b(2); 
        b.setZero();

        double x0 = mesh(elem);
        double x1 = mesh(elem+1);
        double h = x1 - x0;

        // Loop over quadrature points
        for(int q=0; q < nQuadraturePoints; ++q) {
            double qPoint = quad_points[q];
            double w = quad_weights[q];

            // Map reference [0,1] to physical element
            double x = x0 + h * qPoint;

            // Shape functions and gradients at quadrature point
            double phi[2] = {1.0 - qPoint, qPoint};
            double dphi[2] = {-1.0/h, 1.0/h}; // d/dx of shape functions

            // Evaluate coefficients at x
            double diff = diffusion_term(x);
            double transp = transport_term(x);
            double react = reaction_term(x);
            double forc = forcing_term(x);

            for(int i=0; i<2; ++i) {
                b(i) += forc * phi[i] * w * h;
                for(int j=0; j<2; ++j) {
                    mat(i,j) += (
                        diff * dphi[i] * dphi[j] +
                        transp * dphi[j] * phi[i] +
                        react * phi[i] * phi[j]
                    ) * w * h;
                }
            }
        }

        rhs[elem] += b(0);
        rhs[elem+1] += b(1);

        for(int r=0; r<2; r++)
            for(int c=0; c<2; c++)
                triplets.emplace_back(elem+r, elem+c, mat(r, c));
    }
    A.setFromTriplets(triplets.begin(), triplets.end());


    //TODO
    // da spostare fuori
    //TODO
    //capire se si può fare diverso da Point<1>(0.0)
    std::vector<BoundaryCondition<1,1>> boundary_conditions = {BoundaryCondition<1,1>(0, BCType::DIRICHLET, Function<1,1>([](Point<1> p) { return Point<1>(0.0); })),
                                                                  BoundaryCondition<1,1>(1, BCType::NEUMANN, Function<1,1>([](Point<1> p) { return Point<1>(0.0); }))};
    BoundaryConditions<1,1> boundary_conds(boundary_conditions);
    boundary_conds.apply(mesh, A, rhs);

    // // Dirichlet/Neumann boundary handling (unchanged)
    // if (!boundary_conds[0].isNeumann()) {
    //     A.coeffRef(0,0) = 1;
    //     A.coeffRef(0,1) = 0;
    //     rhs[0] = boundary_conds[0].getBoundary()(mesh(0));
    // } else {
    //     rhs[0] -= boundary_conds[0].getBoundary()(mesh(0)) * 0.5;
    // }
    // int n = A.rows()-1;
    // if (!boundary_conds[1].isNeumann()) {
    //     A.coeffRef(n,n-1) = 0;
    //     A.coeffRef(n,n) = 1;
    //     rhs[n] = boundary_conds[1].getBoundary()(mesh.getEnd());
    // } else {
    //     rhs[n] += boundary_conds[1].getBoundary()(mesh.getEnd()) * 0.5;
    // }

    std::cout << "Vector RHS:\n" << rhs << std::endl;
    std::cout << "Matrix A (LHS):\n" << A << std::endl;
}

const char* Fem1D::solverInfoToString(Eigen::ComputationInfo info){
    switch(info){
        case Eigen::ComputationInfo::Success:
            return "Success";
        case Eigen::ComputationInfo::NumericalIssue:
            return "Numerical Issue (prerequisites not satisfied by data)";
        case Eigen::ComputationInfo::NoConvergence:
            return "Did not Converge, max iterations reached";
        case Eigen::ComputationInfo::InvalidInput:
            return "Input invalid";
        default:
            return "Unknown";
    }
}

void Fem1D::solve() {
    Thomas solver;
    try {
        // throw std::runtime_error("Testing BiCG");
        sol = solver.solve(A, rhs);
    } catch (const std::runtime_error& e) {
        std::cout << "Caught exception: " << e.what() << " Solving with BiCGSTAB" << '\n';
        //solve with BiCGSTAB instead
        Eigen::BiCGSTAB<SparseMat> solver;
        solver.setMaxIterations(1e4);
        solver.setTolerance(1e-9);
        solver.compute(A);
        sol = solver.solve(rhs);

        std::cout << "BiCG Status: " << solverInfoToString(solver.info()) << std::endl;
        std::cout << "Iterations: " << solver.iterations() << std::endl;
        std::cout << "Error: " << solver.error() << std::endl;
    }
};

void Fem1D::solve(std::ofstream &fout) {
    solve();
    fout << "x,f(x)" << std::endl;
    for(int i=0 ; i<mesh.getN() ; i++) {
        fout << mesh(i) << "," << sol[i] << std::endl;
    }
};