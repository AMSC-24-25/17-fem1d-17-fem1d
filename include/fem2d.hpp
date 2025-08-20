#ifndef FEM2D_HPP
#define FEM2D_HPP

#include "function.hpp"
#include "grid2D.hpp"
#include "boundary_conditions.hpp"
#include "quadrature.hpp"

#include <fstream>
#include <iostream>
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

// Struttura per boundary conditions moderne
template<unsigned int dim>
class Fem {
private:
    Grid<dim> mesh;
    Function<dim, 1> forcing_term;
    Function<dim, 1> reaction_term;
    Function<dim, 1> diffusion_term;
    Function<dim, 1> transport_term;

    BoundaryConditions<dim, 1> boundaryConditions;

    SparseMat A;
    VectorXd rhs;
    VectorXd solution;

public:
    // Costruttore moderno con BoundaryConditions
    Fem(Grid<dim> grid, Function<dim, 1> forcing, Function<dim, 1> diffusion, 
        Function<dim, 1> transport, Function<dim, 1> reaction,
        const BoundaryConditions<dim, 1>& boundaryConditions);

    // Assemblaggio matrici
    void assemble();
    
    // Risoluzione sistema
    void solve();
    
    // Getter per la soluzione
    const VectorXd& getSolution() const { return solution; }

    void outputVtu(const std::string& filename) const;
    void outputCsv(const std::string& filename) const;
private:
    OrderTwoQuadrature getQuadratureRule() const;

    // Metodi helper per assemblaggio
    void assembleElement(int elemIndex, BarycentricQuadRule& quad, 
                        std::vector<Triplet>& triplets);
        
    void applyDirichletBC();
};

#endif // FEM2D_HPP
