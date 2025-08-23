/// @file boundary_conditions.tpp
/// @brief Template implementation for BoundaryConditions_td classes
/// @details This file contains the implementation of boundary conditions
///          for 1D, 2D and 3D finite element problems

#include "boundary_conditions_td.hpp"
#include "quadrature.hpp"
#include <iostream>

// =============================================================================
// GENERIC TEMPLATE IMPLEMENTATION (2D/3D)
// =============================================================================

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
BoundaryConditions_td<dim, returnDim>::BoundaryConditions_td(
    const std::vector<BoundaryCondition_td<dim, returnDim>>& conditions) 
    : conditions(conditions) {}

// -----------------------------------------------------------------------------
// Apply boundary conditions (main method)
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions_td<dim, returnDim>::apply(const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) {
    
    if (conditions.empty()) {
        std::cout << "No boundary condition to apply" << std::endl;
        return;
    }

    std::cout << "Applying " << conditions.size() << " boundary conditions..." << std::endl;

    std::vector<BoundaryCondition_td<dim, returnDim>> dirichletConditions;

    for (const BoundaryCondition_td<dim, returnDim>& condition : conditions) {
        switch (condition.getType()) {
            case BCType::DIRICHLET:
                dirichletConditions.push_back(condition);
                break;
            case BCType::NEUMANN:
                applyNeumann(condition, mesh, A, rhs, t);
                break;
        }
    }

    for (const BoundaryCondition_td<dim, returnDim>& dirichletCondition : dirichletConditions) {
        applyDirichlet(dirichletCondition, mesh, A, rhs, t);
    }

    std::cout << "Boundary conditions applied successfully" << std::endl;
}

// -----------------------------------------------------------------------------
// Helper methods for applying specific boundary condition types
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions_td<dim, returnDim>::applyDirichlet(
    const BoundaryCondition_td<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) {
    
    IndexVector boundaryNodes = mesh.getBoundaryNodesByTag(bc.getBoundaryId());
    std::cout << "  Applying Dirichlet condition on tag " << bc.getBoundaryId() 
              << " (" << boundaryNodes.size() << " nodes)" << std::endl;

    // This loop is embarassingly parallel: boundaryNodes does not contain duplicates
    #pragma omp parallel for
    for (unsigned int nodeIndex : boundaryNodes) {
        const Point<dim>& nodePoint = mesh.getNode(nodeIndex);
        
        // Set row to zero (only iterates over nonzero elements in the row)
        for (SparseMat::InnerIterator it(A, nodeIndex); it; ++it) {
            it.valueRef() = 0.0;
        }
        
        // Set 1 on diagonal
        A.coeffRef(nodeIndex, nodeIndex) = 1.0;
        // Set dirichlet value in rhs
        rhs[nodeIndex] = bc.getBoundaryFunction(t).value(nodePoint);
    }
}
