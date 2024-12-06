#include "../include/fem1d.hpp"
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

void Fem1D::assemble() {
    FunctionVector phiVect = mesh.getPhiFunctions();
    std::vector<Triplet> triplets;

    triplets.reserve(3 * mesh.getN() - 2);
    
    for(int i=0; i<mesh.getN()-1 ; i++) {
        MatrixXd mat(2, 2);
        VectorXd b(2);
        for(int k1=0 ; k1<=1 ; k1++) {
            for(int k2=0 ; k2<=1 ; k2++) {
                TwoPointsQuadrature quad(
                    (diffusion_term * phiVect[i+k1].getGrad() * phiVect[i+k2].getGrad()) +
                    (reaction_term * phiVect[i+k1] * phiVect[i+k2])
                );

                mat(k1, k2) = quad.integrate(mesh(i), mesh(i+1));
            }
            TwoPointsQuadrature quad2(forcing_term *phiVect[i]);
            b[k1]=quad2.integrate(mesh(i+k1), mesh(i+1+k1));
        }
        // std::cout<< "Small matrix coming from " << i << "th element:\n" << mat << std::endl;
        
        rhs[i] += b[0];
        rhs[i+1] += b[1];
        for (BoundaryCond boundary_cond : boundary_conds){
            if ((boundary_cond).isNeumann()) {
                rhs[i] +=  boundary_cond.getBoundary()(mesh.getEnd()) * phiVect[i](mesh.getEnd());
            }
        }
        
        triplets.emplace_back(i, i, mat(0, 0));
        triplets.emplace_back(i, i + 1, mat(0, 1));
        triplets.emplace_back(i + 1, i, mat(1, 0));
        triplets.emplace_back(i + 1, i + 1, mat(1, 1));

    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Dirichlet
    if (!boundary_conds[0].isNeumann()) {
        A.coeffRef(0,0) = 1;
        A.coeffRef(0,1) = 0;
        rhs[0] = boundary_conds[0].getBoundary()(mesh(0));
    }
    if (!boundary_conds[1].isNeumann()) {
        int n = A.rows()-1;
        A.coeffRef(n,n-1) = 0;
        A.coeffRef(n,n) = 1;
        rhs[n] = boundary_conds[1].getBoundary()(mesh.getEnd());
    }

    // std::cout << "Vector RHS:\n" << rhs << std::endl;
    // std::cout << "Matrix A (LHS):\n" << A << std::endl;
};

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
        sol = solver.ThomasAlgorithm(A, rhs);
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