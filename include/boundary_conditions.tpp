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
        std::cout << "Nessuna condizione al contorno da applicare" << std::endl;
        return;
    }
    
    std::cout << "Applicazione " << conditions.size() << " condizioni al contorno..." << std::endl;

    for (const BoundaryCondition<dim, returnDim>& condition : conditions) {
        switch (condition.getType()) {
            case BCType::DIRICHLET:
                applyDirichlet(condition, mesh, A, rhs);
                break;
            case BCType::NEUMANN:
                applyNeumann(condition, mesh, A, rhs);
                break;
        }
    }
    
    std::cout << "Condizioni al contorno applicate con successo" << std::endl;
}

// -----------------------------------------------------------------------------
// Helper methods for applying specific boundary condition types
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyDirichlet(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    IndexVector boundaryNodes = mesh.getBoundaryNodesByTag(bc.getBoundaryId());
    std::cout << "  Applicando condizione di Dirichlet su tag " << bc.getBoundaryId() 
              << " (" << boundaryNodes.size() << " nodi)" << std::endl;
    
    for (unsigned int nodeIndex : boundaryNodes) {
        const Point<dim>& nodePoint = mesh.getNode(nodeIndex);
        
        // Set row to zero
        for (int k = 0; k < A.outerSize(); ++k) {
            for (SparseMat::InnerIterator it(A, k); it; ++it) {
                if (it.row() == nodeIndex)
                    it.valueRef() = 0.0;
            }
        }
        
        // Set 1 on diagonal
        A.coeffRef(nodeIndex, nodeIndex) = 1.0;
        // Set dirichlet value in rhs
        rhs[nodeIndex] = bc.getBoundaryFunction().value(nodePoint);
    }
}

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyNeumann(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    // Ottieni gli edge di bordo con il tag fisico specifico
    auto boundaryEdges = mesh.getBoundaryEdgesByTag(bc.getBoundaryId());
    std::cout << "  Applicando condizione di Neumann su tag " << bc.getBoundaryId()
              << " (" << boundaryEdges.size() << " edge)" << std::endl;
    
    // Inizializza la quadratura 1D
    GaussLegendre1D quadrature;
    
    // Itera su tutti gli edge di bordo con questo tag
    for (const auto& edge : boundaryEdges) {
        // Calcola i contributi ai nodi dell'edge usando la quadratura
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(edge, bc.getBoundaryFunction(), contributions);
        
        // Aggiungi i contributi al vettore RHS
        const auto& nodeIndices = edge.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
        
        std::cout << "    Edge con nodi [" << nodeIndices[0] << ", " << nodeIndices[1] 
                 << "] - contributi: [" << contributions[0] << ", " << contributions[1] << "]" << std::endl;
    }
}

// =============================================================================
// 1D SPECIALIZATION IMPLEMENTATION
// =============================================================================

// -----------------------------------------------------------------------------
// Apply boundary conditions (1D main method)
// -----------------------------------------------------------------------------

inline void BoundaryConditions<1,1>::apply(const Grid1D& mesh,
    SparseMat& A, VectorXd& rhs) {
    
    if (conditions.empty()) {
        std::cout << "Nessuna condizione al contorno 1D" << std::endl;
        return;
    }
    
    std::cout << "Applicazione " << conditions.size() << " condizioni al contorno 1D..." << std::endl;
    
    for (const auto& bc : conditions) {
        if (bc.getType() == BCType::DIRICHLET) {
            applyDirichlet(bc, mesh, A, rhs);
        } else { // NEUMANN
            applyNeumann(bc, mesh, A, rhs);
        }
    }
    
    std::cout << "Condizioni al contorno 1D applicate con successo" << std::endl;
}

// -----------------------------------------------------------------------------
// Helper methods for applying specific boundary condition types (1D)
// -----------------------------------------------------------------------------

inline void BoundaryConditions<1, 1>::applyDirichlet(
    const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    std::cout << "  Applicando condizione di Dirichlet 1D su tag " << bc.getBoundaryId() << std::endl;
    
    if (bc.getBoundaryId() == 0) {
        // Boundary condition at the left endpoint (x = 0)
        int startIndex = mesh.getStart();
        A.coeffRef(startIndex, startIndex) = 1.0;
        if (startIndex + 1 < A.cols()) {
            A.coeffRef(startIndex, startIndex + 1) = 0.0;
        }
        rhs[startIndex] = bc.getBoundaryFunction().value(Point<1>(mesh.getStart()));
        
    } else if (bc.getBoundaryId() == 1) {
        // Boundary condition at the right endpoint (x = L)
        int endIndex = mesh.getEnd();
        A.coeffRef(endIndex, endIndex) = 1.0;
        if (endIndex - 1 >= 0) {
            A.coeffRef(endIndex, endIndex - 1) = 0.0;
        }
        rhs[endIndex] = bc.getBoundaryFunction().value(Point<1>(mesh.getEnd()));
        
    } else {
        std::cerr << "ERRORE: Tag fisico non valido per condizione di Dirichlet 1D: " 
                  << bc.getBoundaryId() << " (attesi: 0 o 1)" << std::endl;
    }
}

inline void BoundaryConditions<1,1>::applyNeumann(
    const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    std::cout << "  Applicando condizione di Neumann 1D su tag " << bc.getBoundaryId() << std::endl;
    
    if (bc.getBoundaryId() == 0) {
        // Neumann condition at the left endpoint (x = 0)
        // ∂u/∂n = -∂u/∂x at x = 0 (outward normal points left)
        int startIndex = mesh.getStart();
        double neumannValue = bc.getBoundaryFunction().value(Point<1>(mesh.getStart()));
        rhs[startIndex] -= neumannValue * 0.5; // Factor 0.5 comes from FEM integration
        
    } else if (bc.getBoundaryId() == 1) {
        // Neumann condition at the right endpoint (x = L)
        // ∂u/∂n = ∂u/∂x at x = L (outward normal points right)
        int endIndex = mesh.getEnd();
        double neumannValue = bc.getBoundaryFunction().value(Point<1>(mesh.getEnd()));
        rhs[endIndex] += neumannValue * 0.5; // Factor 0.5 comes from FEM integration
        
    } else {
        std::cerr << "ERRORE: Tag fisico non valido per condizione di Neumann 1D: " 
                  << bc.getBoundaryId() << " (attesi: 0 o 1)" << std::endl;
    }
}

// =============================================================================
// 3D TODO SECTION
// =============================================================================

/*
 * TODO: Implementazione per 3D
 * 
 * Per supportare completamente il 3D, è necessario implementare:
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
 * Il pattern sarà simile al 2D:
 * 
 * void applyNeumann3D(const BoundaryCondition<3,1>& bc, const Grid3D& mesh, 
 *                     SparseMat& A, VectorXd& rhs) {
 *     auto boundaryFaces = mesh.getBoundaryFacesByTag(bc.getBoundaryId());
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
