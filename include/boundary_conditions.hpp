#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "function.hpp"
#include "point.hpp"
#include "grid1D.hpp"
#include "grid.hpp"
#include <Eigen/Sparse>
#include <vector>

using Eigen::VectorXd;
using SparseMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// =============================================================================
// ENUMS AND BASIC STRUCTURES
// =============================================================================

/// Enum per il tipo di condizione al contorno
enum class BCType {
    DIRICHLET,  ///< u = valore specificato
    NEUMANN     ///< du/dn = valore specificato
};

// Single boundary condition
template<unsigned int dim, unsigned int returnDim>
class BoundaryCondition {
private:
    int boundaryId;                    // Boundary id
    BCType type;                       // Condition type (Dirichlet/Neumann)
    Function<dim, returnDim> boundaryFunction; // Function defining the condition

public:
    // Constructor
    BoundaryCondition(int tag, BCType bcType, Function<dim, returnDim> func) 
        : boundaryId(tag), type(bcType), boundaryFunction(func) {}
    
    // Constructor for constant value (for convenience)
    BoundaryCondition(int tag, BCType bcType, double value) 
        : boundaryId(tag), type(bcType), 
          boundaryFunction([value](Point<dim> p) { return value; }) {}

    // Getters
    int getBoundaryId() const { return boundaryId; }
    BCType getType() const { return type; }
    Function<dim, returnDim> getBoundaryFunction() const { return boundaryFunction; }
};

// Manages a set of boundary conditions
template <unsigned int dim, unsigned int returnDim>
class BoundaryConditions {
    using IndexVector = std::vector<unsigned int>;

public:
    // Costruttori
    BoundaryConditions() = default;
    BoundaryConditions(const std::vector<BoundaryCondition<dim, returnDim>>& conditions);

    // Metodi per aggiungere condizioni
    void addDirichlet(int boundaryId, Function<dim, returnDim> func) {
        conditions.emplace_back(boundaryId, BCType::DIRICHLET, func);
    }
    void addDirichlet(int boundaryId, Point<returnDim> value) {
        conditions.emplace_back(boundaryId, BCType::DIRICHLET, value);
    }
    void addNeumann(int boundaryId, Function<dim, returnDim> func) {
        conditions.emplace_back(boundaryId, BCType::NEUMANN, func);
    }
    void addNeumann(int boundaryId, Point<returnDim> value) {
        conditions.emplace_back(boundaryId, BCType::NEUMANN, value);
    }

    // Applicazione delle condizioni al contorno
    void apply(const Grid<dim>& mesh, SparseMat& A, VectorXd& rhs);

private:
    std::vector<BoundaryCondition<dim, returnDim>> conditions;

    // Metodi helper per applicazione
    void applyDirichlet(const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
                       SparseMat& A, VectorXd& rhs);
    void applyNeumann(const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
                     SparseMat& A, VectorXd& rhs);
};

// =============================================================================
// TEMPLATE IMPLEMENTATION
// =============================================================================

#include "boundary_conditions.tpp"

#endif // BOUNDARY_CONDITIONS_HPP
