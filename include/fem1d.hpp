#ifndef FEM1D_HPP
#define FEM1D_HPP

#include "function.hpp"
#include "grid1D.hpp"
#include "boundary_cond.hpp"
#include "fem1d.hpp"
#include "quadrature.hpp"
#include "thomas.hpp"

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;
typedef Eigen::Triplet<double> T;

class Fem1D {

    private:
    Grid1D mesh;
    Function forcing_term;
    Function reaction_term;
    Function diffusion_term;
    
    std::vector<BoundaryCond> boundary_conds;

    SparseMat A;
    Eigen::VectorXd rhs;
    Eigen::VectorXd sol;

    public:
    Fem1D(Grid1D mesh, Function forcing_term, Function reaction_term, Function diffusion_term,
        bool isNeuman1, bool isNeuman2, Function value1, Function value2):
        mesh(mesh),
        forcing_term(forcing_term),
        reaction_term(reaction_term),
        diffusion_term(diffusion_term),
        A(mesh.getN(), mesh.getN()),
        rhs(mesh.getN()),
        sol(mesh.getN()) {
            BoundaryCond boundary1(isNeuman1, value1);
            BoundaryCond boundary2(isNeuman2, value2);
            boundary_conds.push_back(boundary1);
            boundary_conds.push_back(boundary2); 
        };
    
    void assemble();
    void solve();
    void solve(std::ofstream &fout);

    inline Eigen::VectorXd getSolution() const {
        return sol;
    };
    
};

#endif  // FEM1D_HPP