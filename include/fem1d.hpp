#ifndef FEM1D_HPP
#define FEM1D_HPP

#include "function.hpp"
#include "grid1D.hpp"
#include "boundary_cond.hpp"
#include "matrix.hpp"

class Fem1D {

    private:
    Grid1D mesh;
    Function forcing_term;
    Function reaction_term;
    Function diffusion_term;
    
    BoundaryConds boundary_conds;

    Matrix A;
    std::vector<double> rhs;

    public:

    Fem1D(int L, int N, Function forcing_term, Function reaction_term, Function diffusion_term,
        bool bc1, bool bc2, Function value1, Function value2):
        mesh(0, L, N),
        forcing_term(forcing_term),
        reaction_term(reaction_term),
        diffusion_term(diffusion_term),
        boundary_conds(bc1, bc2, value1, value2),
        A(N), rhs(N) {};

    void assemble();
    void solve();
    
};



#endif  // FEM1D_HPP