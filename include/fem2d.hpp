#ifndef FEM2D_HPP
#define FEM2D_HPP

#include "function.hpp"
#include "grid2D.hpp"
#include "boundary_cond.hpp"
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
    Function<2> forcing_term;
    Function<2> reaction_term;
    Function<2> diffusion_term;
    Function<2> transport_term;
    
    // Condizioni al contorno (legacy per compatibilità)
    bool isNeumann1, isNeumann2;
    Function<2> boundary1, boundary2;
    
    // Boundary conditions moderne
    BoundaryConditions boundaryConditions;
    
    SparseMat A;
    VectorXd rhs;
    VectorXd solution;

public:
    // Costruttore moderno con BoundaryConditions
    Fem2D(Grid2D grid, Function<2> forcing, Function<2> diffusion, 
          Function<2> transport, Function<2> reaction,
          const BoundaryConditions& boundaryConditions);
    
    // Costruttore legacy (per compatibilità)
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

    void outputVtk(const std::string& filename) const;
private:
    // Metodi helper per assemblaggio
    void assembleElement(int elemIndex, orderTwoQuadrature& quad, 
                        std::vector<Triplet>& triplets);
        
    void applyDirichletBC();
};

#endif // FEM2D_HPP
