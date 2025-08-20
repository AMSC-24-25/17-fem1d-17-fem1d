#ifndef FEM1D_HPP
#define FEM1D_HPP

#include "function.hpp"
#include "grid1D.hpp"
#include "boundary_conditions.hpp"
#include "thomas.hpp"

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>

using Eigen::VectorXd;
using Eigen::MatrixXd;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;
typedef Eigen::Triplet<double> Triplet;

class Fem1D {

    private:
    Grid1D mesh;
    Function<1, 1> forcing_term;
    Function<1, 1> reaction_term;
    Function<1, 1> diffusion_term;
    Function<1, 1> transport_term;

    //TODO
    //porta qui le boundary e passale
    BoundaryConditions<1, 1> boundaryConditions;


    SparseMat A;
    VectorXd rhs;
    VectorXd sol;

    public:
    Fem1D(Grid1D mesh, Function<1, 1> forcing_term, Function<1, 1> diffusion_term, 
            Function<1, 1> transport_term, Function<1, 1> reaction_term, 
            BoundaryConditions<1, 1> boundaryConditions):
        mesh(mesh),
        forcing_term(forcing_term),
        diffusion_term(diffusion_term),
        transport_term(transport_term),
        reaction_term(reaction_term),
        A(mesh.getN(), mesh.getN()),
        rhs(mesh.getN()),
        sol(mesh.getN()),
        boundaryConditions(boundaryConditions)
    {};

    void assemble();
    void solve();
    void solve(std::ofstream &fout);

    inline VectorXd getSolution() const {
        return sol;
    };

    private:
    static const char* solverInfoToString(Eigen::ComputationInfo info);
};

#endif  // FEM1D_HPP