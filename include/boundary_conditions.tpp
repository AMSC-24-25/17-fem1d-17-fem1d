/// @file boundary_conditions.tpp
/// @brief Template implementation for BoundaryConditions classes
/// @details This file contains the implementation of boundary conditions
///          for 1D, 2D and 3D finite element problems

#include "boundary_conditions.hpp"
#include "quadrature.hpp"
#include <iostream>

// =============================================================================
// GENERIC TEMPLATE IMPLEMENTATION (2D/3D)
// =============================================================================

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
BoundaryConditions<dim, returnDim>::BoundaryConditions(
    const std::vector<BoundaryCondition<dim, returnDim>>& conditions) 
    : conditions(conditions) {}

// -----------------------------------------------------------------------------
// Apply boundary conditions (main method)
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::apply(const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    if (conditions.empty()) {
        std::cout << "No boundary condition to apply" << std::endl;
        return;
    }

    std::cout << "Applying " << conditions.size() << " boundary conditions..." << std::endl;

    std::vector<BoundaryCondition<dim, returnDim>> dirichletConditions;

    for (const BoundaryCondition<dim, returnDim>& condition : conditions) {
        switch (condition.getType()) {
            case BCType::DIRICHLET:
                dirichletConditions.push_back(condition);
                break;
            case BCType::NEUMANN:
                applyNeumann(condition, mesh, A, rhs);
                break;
        }
    }

    for (const BoundaryCondition<dim, returnDim>& dirichletCondition : dirichletConditions) {
        applyDirichlet(dirichletCondition, mesh, A, rhs);
    }

    std::cout << "Boundary conditions applied successfully" << std::endl;
}

// -----------------------------------------------------------------------------
// Helper methods for applying specific boundary condition types
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyDirichlet(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
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
        rhs[nodeIndex] = bc.getBoundaryFunction().value(nodePoint);
    }
}

// =============================================================================
// 3D TODO SECTION
// =============================================================================

/*
 * TODO: 3D implementation
 * 
 * To fully support 3D, you need to implement:
 * 
 * 1. Grid3D class with tetrahedra/hexahedra
 * 2. BoundaryCell<2> for triangular/quadrilateral boundary faces
 * 3. GaussLegendre2D quadrature for surface integration
 * 4. Geometric functions:
 *    - mapToGlobalFace()
 *    - getShapeFunctions2D()
 *    - faceArea()
 * 5. Specialization BoundaryConditions<3,1>
 * 
 * The pattern will be similar to 2D:
 * 
 * void applyNeumann3D(const BoundaryCondition<3,1>& bc, const Grid3D& mesh, 
 *                     SparseMat& A, VectorXd& rhs) {
 *     auto boundaryFaces = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
 *     GaussLegendre2D quadrature;
 *     
 *     for (const auto& face : boundaryFaces) {
 *         std::vector<double> contributions;
 *         quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
 *         
 *         const auto& nodeIndices = face.getNodeIndexes();
 *         for (size_t i = 0; i < nodeIndices.size(); ++i) {
 *             rhs[nodeIndices[i]] += contributions[i];
 *         }
 *     }
 * }
 */
