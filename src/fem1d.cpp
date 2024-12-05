#include "fem1d.hpp"
#include "quadrature.hpp"
#include "thomas.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

typedef Triplet<double> Triplet;

void Fem1D::assemble() {
    FunctionVector phiVect = mesh.getPhiFunctions();
    std::vector<Eigen::Triplet<double>> triplets;

    triplets.reserve(3 * mesh.getN() - 2);
    
    for(int i=0; i<mesh.getN()-1 ; i++) {

        MatrixXd mat(2, 2);
        VectorXd b(2);
        for(int k1=0 ; k1<=1 ; k1++) {
            for(int k2=0 ; k2<=1 ; k2++) {

                TwoPointsQuadrature quad(
                    diffusion_term *
                    phiVect[i+k1].getGrad() *
                    phiVect[i+k2].getGrad()
                );

                mat(k1, k2) = quad.integrate(mesh(i), mesh(i+1));
            }

        }

        TwoPointsQuadrature quad(forcing_term*phiVect[i]);


        for (BoundaryCond boundary_cond : boundary_conds){
            if ((boundary_cond).isNeumann()) {
                b[i] = quad.integrate(mesh(i), mesh(i+1)) +
                       boundary_cond.getBoundary()(mesh.getEnd()) * phiVect[i](mesh.getEnd());
            }
        }

        TwoPointsQuadrature quad2(forcing_term *phiVect[i]);
        b[i] = quad2.integrate(mesh(i), mesh(i+1));
        
        triplets.emplace_back(i, i, mat(0, 0));
        triplets.emplace_back(i, i + 1, mat(0, 1));
        triplets.emplace_back(i + 1, i, mat(1, 0));
        triplets.emplace_back(i + 1, i + 1, mat(1, 1));

    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Dirichlet
    if (!boundary_conds[0].isNeumann()) {
        A.coeffRef(0,0) = 0;
        A.coeffRef(0,1) = 0;
    }
    if (!boundary_conds[1].isNeumann()) {
        int n = A.rows()-1;
        A.coeffRef(n,n-1) = 0;
        A.coeffRef(n,n) = 0;
    }

};

void Fem1D::solve() {
    Thomas solver;
    VectorXd solution = solver.ThomasAlgorithm(A, rhs);
    for (int i=0; i < solution.size(); i++){
        std::cout << solution[i] << std::endl;
    }
};