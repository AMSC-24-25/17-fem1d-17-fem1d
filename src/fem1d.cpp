#include "fem1d.hpp"


void Fem1D::assemble() {
    FunctionVector phiVect = mesh.getPhiFunctions();
    std::vector<Triplet> triplets;

    std::cout << typeid(phiVect).name() << std::endl;

    triplets.reserve(3 * mesh.getN() - 2);
    
    for(int i=0; i<mesh.getN()-1 ; i++) {
        MatrixXd mat(2, 2);
        VectorXd b(2);
        for(int k1=0 ; k1<=1 ; k1++) {
            for(int k2=0 ; k2<=1 ; k2++) {

                // JxW = Weight of node
                const double JxW = 0.5;
                TwoPointsQuadrature quad(
                    (diffusion_term * phiVect[i+k1].getGrad()[0] * phiVect[i+k2].getGrad()[0] * JxW) +
                    (transport_term * phiVect[i+k1].getGrad()[0] * phiVect[i+k2].getGrad()[0] * JxW) +
                    (reaction_term * phiVect[i+k1].getGrad()[0] * phiVect[i+k2].getGrad()[0] * JxW)
                );

                mat(k1, k2) = quad.integrate(mesh(i), mesh(i+1));
            }
            SimpsonQuadrature quad2(forcing_term * phiVect[i]);
            b[k1]=quad2.integrate(mesh(i+k1), mesh(i+1+k1));
        }
        
        rhs[i] += b[0];
        rhs[i+1] += b[1];

        for(int r=0; r<2; r++)
            for(int c=0; c<2; c++)
                triplets.emplace_back(i+r, i+c, mat(r, c));
    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Dirichlet
    if (!boundary_conds[0].isNeumann()) {
        A.coeffRef(0,0) = 1;
        A.coeffRef(0,1) = 0;
        rhs[0] = boundary_conds[0].getBoundary()(mesh(0));
    } 
    else { // Neumann
        rhs[0] -= boundary_conds[0].getBoundary()(mesh(0)) * 0.5;
    }

    int n = A.rows()-1;
    // Dirichlet
    if (!boundary_conds[1].isNeumann()) {
        A.coeffRef(n,n-1) = 0;
        A.coeffRef(n,n) = 1;
        rhs[n] = boundary_conds[1].getBoundary()(mesh.getEnd());
    } 
    else { // Neumann
        rhs[n] += boundary_conds[1].getBoundary()(mesh.getEnd()) * 0.5;
    }

    std::cout << "Vector RHS:\n" << rhs << std::endl;
    std::cout << "Matrix A (LHS):\n" << A << std::endl;
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