#include "../include/boundary_conditions.hpp"
#include "../include/quadrature.hpp"
#include <iostream>

template <>
void BoundaryConditions<3, 1>::applyNeumann(
    const BoundaryCondition<3, 1>& bc, const Grid<3>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    // Ottieni le facce di bordo con il tag fisico specifico
    //from auto to std::vector<BoundaryCell<2>>
    std::vector<BoundaryCell<2>> boundaryFaces = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applicando condizione di Neumann 3D su tag " << bc.getBoundaryId()
              << " (" << boundaryFaces.size() << " facce)" << std::endl;
    
    // Inizializza la quadratura 2D per facce triangolari
    GaussLegendre2D quadrature;
    
    // Itera su tutte le facce di bordo con questo tag
    //from const auto& to const BoundaryCell<2>&
    for (const BoundaryCell<2>& face : boundaryFaces) {
        // Calcola i contributi ai nodi della faccia usando la quadratura
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
        
        // Aggiungi i contributi al vettore RHS
        //from const auto& to const std::vector<unsigned int>&
        const std::vector<unsigned int>& nodeIndices = face.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
    }
}

// template <unsigned int dim, unsigned int returnDim>
template <>
void BoundaryConditions<2, 1>::applyNeumann(
    const BoundaryCondition<2, 1>& bc, const Grid<2>& mesh, 
    SparseMat& A, VectorXd& rhs) {
    
    // Ottieni gli edge di bordo con il tag fisico specifico
    //from auto to std::vector<BoundaryCell<1>>
    std::vector<BoundaryCell<1>> boundaryEdges = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applicando condizione di Neumann su tag " << bc.getBoundaryId()
              << " (" << boundaryEdges.size() << " edge)" << std::endl;
    
    // Inizializza la quadratura 1D
    GaussLegendre1D quadrature;
    
    // Itera su tutti gli edge di bordo con questo tag
    //from const auto& to const BoundaryCell<1>&
    for (const BoundaryCell<1>& edge : boundaryEdges) {
        // Calcola i contributi ai nodi dell'edge usando la quadratura
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(edge, bc.getBoundaryFunction(), contributions);
        
        // Aggiungi i contributi al vettore RHS
        //from const auto& to const std::vector<unsigned int>&
        const std::vector<unsigned int>& nodeIndices = edge.getNodeIndexes();
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

    const BoundaryCell<0> boundaryCell = mesh.getBoundaryCellsByTag(bc.getBoundaryId())[0];
    unsigned int node_index = boundaryCell.getNodeIndex(0);
    double node = boundaryCell.getNode(0);

    int adjacentNode = (node_index == 0) ? (node_index+1) : (node_index-1);
    double h = mesh.getNode(adjacentNode) - mesh.getNode(node_index);
    h = std::abs(h);

    // ∂u/∂n = -∂u/∂x at x = 0 (outward normal points left)
    // ∂u/∂n = ∂u/∂x at x = L (outward normal points right)
    int sign = (node_index == 0) ? -1 : 1;
    double neumannValue = sign * bc.getBoundaryFunction().value(node);
    rhs[node_index] -= neumannValue * h/2.0;
}
