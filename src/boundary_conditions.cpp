#include "../include/boundary_conditions.hpp"
#include "../include/quadrature.hpp"
#include <iostream>

template <>
void BoundaryConditions<3, 1>::applyNeumann(
    const BoundaryCondition<3, 1>& bc, const Grid<3>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    // Ottieni le facce di bordo con il tag fisico specifico
    auto boundaryFaces = mesh.getBoundaryFacesByTag(bc.getBoundaryId());
    std::cout << "  Applicando condizione di Neumann 3D su tag " << bc.getBoundaryId()
              << " (" << boundaryFaces.size() << " facce)" << std::endl;
    
    // Inizializza la quadratura 2D per facce triangolari
    GaussLegendre2D quadrature;
    
    // Itera su tutte le facce di bordo con questo tag
    for (const auto& face : boundaryFaces) {
        // Calcola i contributi ai nodi della faccia usando la quadratura
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
        
        // Aggiungi i contributi al vettore RHS
        const auto& nodeIndices = face.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
        
        std::cout << "    Faccia con nodi [" << nodeIndices[0] << ", " << nodeIndices[1] 
                  << ", " << nodeIndices[2] << "] - contributi: [" << contributions[0] 
                  << ", " << contributions[1] << ", " << contributions[2] << "]" << std::endl;
    }
}

// template <unsigned int dim, unsigned int returnDim>
template <>
void BoundaryConditions<2, 1>::applyNeumann(
    const BoundaryCondition<2, 1>& bc, const Grid<2>& mesh, 
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
template<>
void BoundaryConditions<1,1>::applyNeumann(
    const BoundaryCondition<1, 1>& bc, const Grid<1>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    std::cout << "  Applicando condizione di Neumann 1D su tag " << bc.getBoundaryId() << std::endl;
    
    unsigned int node_index = mesh.getBoundaryNodesByTag(bc.getBoundaryId())[0];
    double node = mesh.getNode(node_index);

    // ∂u/∂n = -∂u/∂x at x = 0 (outward normal points left)
    // ∂u/∂n = ∂u/∂x at x = L (outward normal points right)
    double sign = (node_index == 0) ? -1.0 : 1.0;
    double neumannValue = sign * bc.getBoundaryFunction().value(node);
    rhs[node_index] -= neumannValue * 0.5; // Factor 0.5 comes from FEM integration
}
