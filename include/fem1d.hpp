#ifndef FEM1D_HPP
#define FEM1D_HPP

#include <vector>
#include <fstream>

#include <iostream>
#include "function.hpp"
#include "grid1D.hpp"
#include "boundary_cond.hpp"
#include "fem1d.hpp"
#include "quadrature.hpp"
#include "thomas.hpp"
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

using Eigen::VectorXd;
using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;
typedef Eigen::Triplet<double> Triplet;

class Fem1D {

    private:
    Grid1D mesh;
    Function forcing_term;
    Function reaction_term;
    Function diffusion_term;
    
    std::vector<BoundaryCond> boundary_conds;

    SparseMat A;
    VectorXd rhs;
    VectorXd sol;

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

    inline VectorXd getSolution() const {
        return sol;
    };
    
};

#endif  // FEM1D_HPP