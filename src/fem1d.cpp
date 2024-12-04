#include "fem1d.hpp"
#include "quadrature.hpp"
#include "../include/matrix.hpp"
#include "thomas.hpp"

void Fem1D::assemble() {

    for(int i=0; i<mesh.getN()-1 ; i++) {

        Matrix mat(2);
        Vector b(2);
        FunctionVector phiVect = mesh.getPhiFunctions();
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
            
        if(boundary_conds.getBc1())
            b[i] = quad.integrate(mesh(i), mesh(i+1)) + 
                    boundary_conds.getValue1()(mesh.getEnd()) * phiVect[i](mesh.getEnd());

        if(boundary_conds.getBc2())
            b[i] = quad.integrate(mesh(i), mesh(i+1)) + 
                    boundary_conds.getValue2()(mesh.getEnd()) * phiVect[i](mesh.getEnd());

        TwoPointsQuadrature quad2(forcing_term *phiVect[i]);
        b[i] = quad2.integrate(mesh(i), mesh(i+1));
        
        A.add_contribution(i, mat);

    }
    

    // Dirichlet
    if (!boundary_conds.getBc1()) {
        A(0,0) = 0;
        A(0,1) = 0;
    }
    if (!boundary_conds.getBc2()) {
        int n = A.getSize();
        A(n,n-1) = 0;
        A(n,n) = 0;
    }

};

void Fem1D::solve() {
    Thomas solver;
    Vector solution = solver.ThomasAlgorithm(A, rhs);
};