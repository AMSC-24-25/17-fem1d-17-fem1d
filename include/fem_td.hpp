/**
 * @file fem_td.hpp
 * @brief Time-dependent finite element solver using Theta-method
 * 
 * Implements finite element solution for time-dependent PDEs:
 * 
 * Uses Theta-method for time discretization with configurable Theta parameter.
 */
#pragma once
#include <vector>
#include <string>
#include <Eigen/Sparse>

#include "grid.hpp"
#include "cell.hpp"
#include "point.hpp"
#include "function.hpp"
#include "quadrature.hpp"
#include "boundary_conditions_td.hpp"

using SparseMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Triplet   = Eigen::Triplet<double>;
using MatrixXd  = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXd  = Eigen::VectorXd;

/**
 * @brief Time-dependent FEM solver with Theta-method time stepping
 */
template<unsigned int dim>
class FemTD {
public:
    using ForcingTD = std::function<double(const Point<dim>&, double)>;

    FemTD(Grid<dim> grid,
          Function<dim,1> diffusion,
          Function<dim,dim> transport,
          Function<dim,1> reaction,
          const BoundaryConditions_td<dim,1>& bc,
          QuadratureRule<dim> quadrature);

    // Configuration methods
    void set_forcing(ForcingTD f_td);                 // Set forcing function f(x,t)
    void set_initial_condition(Function<dim,1> u0);   // Set initial condition u0(x)

    double step(double t_new, double dt, double theta);  // Single time step

    // Complete simulation over time interval [0,T]
    void run(double T, double dt, double theta,
             const std::string& vtu_prefix = "",
             const std::string& csv_prefix = "");

    // Output methods compatible with steady-state Fem<dim>
    void outputCsv(const std::string& filename) const;
    void outputVtu(const std::string& filename) const;

    const VectorXd& solution_vec() const { return u_; }
    const Grid<dim>& mesh() const { return mesh_; }

private:
    Grid<dim>                 mesh_;
    Function<dim,1>           diffusion_;    // μ(x) - diffusion coefficient
    Function<dim,dim>         transport_;    // b(x) - transport coefficient vector
    Function<dim,1>           reaction_;     // r(x) - reaction coefficient
    BoundaryConditions_td<dim,1> bc_;        // Time-dependent boundary conditions
    mutable QuadratureRule<dim>       quad_; // Quadrature rule (mutable for const methods)

    ForcingTD                 forcing_td_;   // f(x,t) - time-dependent forcing
    Function<dim,1>           u0_;           // u0(x) - initial condition

    SparseMat M_, K_, A_;                    // Mass, stiffness, and system matrices
    VectorXd  rhs_;                          // Right-hand side vector
    VectorXd  u_, u_old_;                    // Current and previous solutions

private:
    // Thread-local storage for OpenMP parallelization
#ifdef _OPENMP
    std::vector<VectorXd> F_threads;                              // Per-thread forcing vectors
    std::vector<std::vector<Triplet>> tripletM_thr, tripletK_thr; // Per-thread triplet storage
    int nthreads;                                                 // Number of OpenMP threads
#endif
    std::vector<Triplet> tripletM, tripletK;           // Matrix assembly triplets
    VectorXd f_new, f_old;                             // Forcing vectors at current and previous step

    // Assembly methods
    void assemble_time_invariant();                     // Assemble mass and stiffness matrices
    void assemble_M_and_K_element(int elem,             // Helper method
                                  std::vector<Triplet>& tripM,
                                  std::vector<Triplet>& tripK); 
    void build_load(double t);                          // Build forcing_new vector at time t

    void apply_initial_condition();                     // Apply initial condition u0
};
