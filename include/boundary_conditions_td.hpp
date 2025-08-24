#ifndef BOUNDARY_CONDITIONS_TD_HPP
#define BOUNDARY_CONDITIONS_TD_HPP

#include "function.hpp"
#include "bctype.hpp"
#include "point.hpp"
#include "grid1D.hpp"
#include "grid.hpp"
#include <Eigen/Sparse>
#include <vector>

using Eigen::VectorXd;
using SparseMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// Single boundary condition
template<unsigned int dim, unsigned int returnDim>
class BoundaryCondition_td {
private:
    int boundaryId;                    // Boundary id
    BCType type;                       // Condition type (Dirichlet/Neumann)
    fun_td<dim, returnDim> boundaryFunction; // Function defining the condition

public:
    // Constructor
    BoundaryCondition_td(int tag, BCType bcType, fun_td<dim, returnDim> func) 
        : boundaryId(tag), type(bcType), boundaryFunction(func) {}
    
    // Constructor for constant value (for convenience)
    BoundaryCondition_td(int tag, BCType bcType, double value) 
        : boundaryId(tag), type(bcType), 
          boundaryFunction([value](Point<dim> p, double t) { return value; }) {}

    // Getters
    int getBoundaryId() const { return boundaryId; }
    BCType getType() const { return type; }
    Function<dim, returnDim> getBoundaryFunction(double t) const { 
        return castToSteadyFunction(boundaryFunction, t);
    }
    fun_td<dim, returnDim> getBoundaryFunction() const {
        return boundaryFunction;
    }
};

// Manages a set of boundary conditions
template <unsigned int dim, unsigned int returnDim>
class BoundaryConditions_td {
    using IndexVector = std::vector<unsigned int>;

public:
    // Constructors
    BoundaryConditions_td() = default;
    BoundaryConditions_td(const std::vector<BoundaryCondition_td<dim, returnDim>>& conditions);

    // Methods to add conditions
    void addDirichlet(int boundaryId, fun_td<dim, returnDim> func) {
        conditions.emplace_back(boundaryId, BCType::DIRICHLET, func);
    }
    void addDirichlet(int boundaryId, Point<returnDim> value) {
        conditions.emplace_back(boundaryId, BCType::DIRICHLET, value);
    }
    void addNeumann(int boundaryId, fun_td<dim, returnDim> func) {
        conditions.emplace_back(boundaryId, BCType::NEUMANN, func);
    }
    void addNeumann(int boundaryId, Point<returnDim> value) {
        conditions.emplace_back(boundaryId, BCType::NEUMANN, value);
    }

    // Application of boundary conditions
    void apply(const Grid<dim>& mesh, SparseMat& A, VectorXd& rhs, double t);

private:
    std::vector<BoundaryCondition_td<dim, returnDim>> conditions;

    // Helper methods for application
    void applyDirichlet(const BoundaryCondition_td<dim, returnDim>& bc, const Grid<dim>& mesh, 
                       SparseMat& A, VectorXd& rhs, double t);
    void applyNeumann(const BoundaryCondition_td<dim, returnDim>& bc, const Grid<dim>& mesh, 
                     SparseMat& A, VectorXd& rhs, double t);
};

#include "boundary_conditions_td.tpp"

#endif // BOUNDARY_CONDITIONS_HPP
