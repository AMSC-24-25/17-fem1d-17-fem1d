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
 * FEM time-dependent (θ–method) con le tue strutture dati.
 * PDE: u_t - div(μ∇u) + b·∇u + r u = f(x,t)
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

    // Impostazioni
    void set_forcing(ForcingTD f_td);                 // f(x,t)
    void set_initial_condition(Function<dim,1> u0);   // u0(x)

    // Un singolo passo di tempo (costruisce LHS/RHS e risolve)
    double step(double t_new, double dt, double theta);

    // Esecuzione completa su [0,T] con passo dt
    void run(double T, double dt, double theta,
             const std::string& vtu_prefix = "",
             const std::string& csv_prefix = "");

    // Output compatibili con Fem<dim>
    void outputCsv(const std::string& filename) const;
    void outputVtu(const std::string& filename) const;

    const VectorXd& solution_vec() const { return u_; }
    const Grid<dim>& mesh() const { return mesh_; }

private:
    Grid<dim>                 mesh_;
    Function<dim,1>           diffusion_;    // μ(x)
    Function<dim,dim>         transport_;    // b(x)
    Function<dim,1>           reaction_;     // r(x)
    BoundaryConditions_td<dim,1> bc_;
    mutable QuadratureRule<dim>       quad_; // ← Aggiunto mutable qui

    ForcingTD                 forcing_td_;   // f(x,t)
    Function<dim,1>           u0_;           // u0(x)

    SparseMat M_, K_, A_;
    VectorXd  rhs_;
    VectorXd  u_, u_old_;

private:
    // The following attributes are saved at object level for improved memory efficiency
#ifdef _OPENMP
    std::vector<VectorXd> F_threads;
    std::vector<std::vector<Triplet>> tripletM_thr, tripletK_thr;
#endif
    std::vector<Triplet> tripletM, tripletK;

    void assemble_time_invariant();
    void assemble_M_and_K_element(int elem,
                                  std::vector<Triplet>& tripM,
                                  std::vector<Triplet>& tripK);
    void build_load(VectorXd& F, double t) const; // ∫ f(x,t) φ_i

    void apply_initial_condition();
};
