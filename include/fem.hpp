#ifndef FEM_HPP
#define FEM_HPP

#include "function.hpp"
#include "grid.hpp"
#include "boundary_conditions.hpp"
#include "quadrature.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::Vector3d;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;
typedef Eigen::Triplet<double> Triplet;

// Structure for modern boundary conditions
template<unsigned int dim>
class Fem {
private:
    Grid<dim> mesh;
    Function<dim, 1> forcing_term;
    Function<dim, 1> reaction_term;
    Function<dim, 1> diffusion_term;
    Function<dim, dim> transport_term;

    BoundaryConditions<dim, 1> boundaryConditions;
    QuadratureRule<dim> quadrature;

    SparseMat A;
    VectorXd rhs;
    VectorXd solution;

public:
    // Modern constructor with BoundaryConditions
    Fem(Grid<dim> grid, Function<dim, 1> forcing, Function<dim, 1> diffusion, 
        Function<dim, dim> transport, Function<dim, 1> reaction,
        const BoundaryConditions<dim, 1>& boundaryConditions, QuadratureRule<dim> quadrature);

    // Matrix assembly
    void assemble();
    
    // System solution
    void solve();
    
    // Getter for the solution
    const VectorXd& getSolution() const { return solution; }

    void outputVtu(const std::string& filename) const;
    void outputCsv(const std::string& filename) const;

private:

    // Helper methods for assembly
    void assembleElement(int elemIndex, std::vector<Triplet>& triplets, VectorXd& local_rhs) const;
};

#endif // FEM2D_HPP
