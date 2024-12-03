#ifndef FEM1D_HPP
#define FEM1D_HPP

#include <vector>

#include "function.hpp"
#include "grid1D.hpp"
#include "boundary_cond.hpp"
#include "matrix.hpp"
#include "vector.hpp"

class Fem1D {

    private:
    Grid1D mesh;
    Function forcing_term;
    Function reaction_term;
    Function diffusion_term;
    
    std::vector<BoundaryCond> boundary_conds;

    Matrix A;
    Vector rhs;

    public:
    Fem1D(int L, int N, Function forcing_term, Function reaction_term, Function diffusion_term, bool isNeuman1, bool isNeuman2, Function value1, Function value2):
        mesh(0, L, N),
        forcing_term(forcing_term),
        reaction_term(reaction_term),
        diffusion_term(diffusion_term),
        A(N), rhs(N) {
            BoundaryCond boundary1(isNeuman1, value1);
            BoundaryCond boundary2(isNeuman2, value2);
            boundary_conds.push_back(boundary1);
            boundary_conds.push_back(boundary2);

    };

    void assemble();
    void solve();
    
};



#endif  // FEM1D_HPP