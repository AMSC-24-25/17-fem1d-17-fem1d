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

            }
            //TODO: RHS assemble
            // controlla neumann

            /*if (boundary_conds.getBc1() && i == 0){
                TwoPointsQuadrature quad2(
                        forcing_term *
                        phis[i+k1] +
                        phis[i+k1].value(0) );

                rhs[i] =
            }*/
            TwoPointsQuadrature quad2(forcing_term *phis[i+k1]);
            rhs[i] = quad2.integrate(mesh(i), mesh(i+1));
        }
        A.add_contributions(mat, i);

        // Aggiusta Dirichlet

    }
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
};