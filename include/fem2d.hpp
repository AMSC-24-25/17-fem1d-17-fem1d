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
class Fem2D {
private:
    Grid2D mesh;
    Function<2, 1> forcing_term;
    Function<2, 1> reaction_term;
    Function<2, 1> diffusion_term;
    Function<2, 1> transport_term;

    BoundaryConditions<2, 1> boundaryConditions;

    SparseMat A;
    VectorXd rhs;
    VectorXd solution;

public:
    // Costruttore moderno con BoundaryConditions
    Fem2D(Grid2D grid, Function<2, 1> forcing, Function<2, 1> diffusion, 
          Function<2, 1> transport, Function<2, 1> reaction,
          const BoundaryConditions<2, 1>& boundaryConditions);

    // Assemblaggio matrici
    void assemble();
    
    // Risoluzione sistema
    void solve(std::ofstream& output);
    
    // Getter per la soluzione
    const VectorXd& getSolution() const { return solution; }

    void outputVtk(const std::string& filename) const;
private:
    // Metodi helper per assemblaggio
    void assembleElement(int elemIndex, BarycentricQuadRule& quad, 
                        std::vector<Triplet>& triplets);
        
    void applyDirichletBC();
};

#endif // FEM2D_HPP
