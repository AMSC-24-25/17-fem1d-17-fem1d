#include "fem1d.hpp"
#include "quadrature.hpp"

void Fem1D::assemble() {

    for(int i=0; i<mesh.getN()-1 ; i++) {

        Matrix mat(2);
        FunctionVector phis = mesh.getPhiFunctions();
        for(int k1=0 ; k1<=1 ; k1++) {
            for(int k2=0 ; k2<=1 ; k2++) {

                TwoPointsQuadrature quad(
                    diffusion_term *
                    phis[i+k1].getGrad() *
                    phis[i+k2].getGrad());

                mat(k1, k2) = quad.integrate(mesh(i), mesh(i+1));

                //TODO: RHS assemble
                // controlla neumann
            }
        }

        // Aggiusta Dirichlet

    }

};

void Fem1D::solve() {
};