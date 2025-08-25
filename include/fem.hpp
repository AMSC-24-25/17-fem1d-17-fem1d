/**
 * @file fem.hpp
 * @brief Steady-state finite element solver
 */
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

/**
 * @brief Finite element solver for steady-state PDEs
 */
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
    // Constructor with all PDE coefficients and boundary conditions
    Fem(Grid<dim> grid, Function<dim, 1> forcing, Function<dim, 1> diffusion, 
        Function<dim, dim> transport, Function<dim, 1> reaction,
        const BoundaryConditions<dim, 1>& boundaryConditions, QuadratureRule<dim> quadrature);

    void assemble();  // Assemble system matrix and RHS
    void solve();     // Solve linear system
    
    const VectorXd& getSolution() const { return solution; }

    // Output methods
    void outputVtu(const std::string& filename) const;  // VTU format for ParaView
    void outputCsv(const std::string& filename) const;  // CSV format for analysis

private:
    // Element assembly helper
    void assembleElement(int elemIndex, std::vector<Triplet>& triplets, VectorXd& local_rhs) const;
};

#endif // FEM_HPP
