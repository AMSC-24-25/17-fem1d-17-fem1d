#ifndef FEM2D_HPP
#define FEM2D_HPP

#include "function.hpp"
#include "grid2D.hpp"
#include "boundary_cond.hpp"
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

class Fem2D {
private:
    Grid2D mesh;
    Function<2> forcing_term;
    Function<2> reaction_term;
    Function<2> diffusion_term;
    Function<2> transport_term;
    
    // Condizioni al contorno (semplified per ora)
    bool isNeumann1, isNeumann2;
    Function<2> boundary1, boundary2;
    
    SparseMat A;
    VectorXd rhs;
    VectorXd solution;

public:
    // Costruttore
    Fem2D(Grid2D grid, Function<2> forcing, Function<2> diffusion, 
          Function<2> transport, Function<2> reaction,
          bool isNeumann1, bool isNeumann2, 
          Function<2> boundary1, Function<2> boundary2);
    
    // Assemblaggio matrici
    void assemble();
    
    // Risoluzione sistema
    void solve(std::ofstream& output);
    
    // Getter per la soluzione
    const VectorXd& getSolution() const { return solution; }

private:
    // Metodi helper per assemblaggio
    void assembleElement(int elemIndex, BarycentricQuadRule& quad, 
                        std::vector<Triplet>& triplets);
        
    void applyBoundaryConditions();
};

#endif // FEM2D_HPP
