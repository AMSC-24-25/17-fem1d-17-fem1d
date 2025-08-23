#include "boundary_conditions_td.hpp"
#include "quadrature.hpp"
#include <iostream>

template <>
void BoundaryConditions_td<3, 1>::applyNeumann(
    const BoundaryCondition_td<3, 1>& bc, const Grid<3>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) {
    
    // Get boundary faces with the specified physical tag
    //from auto to std::vector<BoundaryCell<2>>
    std::vector<BoundaryCell<2>> boundaryFaces = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applying Neumann condition on tag " << bc.getBoundaryId()
              << " (" << boundaryFaces.size() << " faces)" << std::endl;

    // Initialize 2D quadrature for triangular faces
    GaussLegendre2D quadrature;
    
    // Iterate over all boundary faces with this tag
    //from const auto& to const BoundaryCell<2>&
    for (const BoundaryCell<2>& face : boundaryFaces) {
    // Compute the contributions to the face nodes using quadrature
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(t), contributions);
        
    // Add the contributions to the RHS vector
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
void BoundaryConditions_td<2, 1>::applyNeumann(
    const BoundaryCondition_td<2, 1>& bc, const Grid<2>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) {
    
    // Get boundary edges with the specified physical tag
    //from auto to std::vector<BoundaryCell<1>>
    std::vector<BoundaryCell<1>> boundaryEdges = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applying Neumann condition on tag " << bc.getBoundaryId()
              << " (" << boundaryEdges.size() << " edge)" << std::endl;
    
    // Initialize 1D quadrature
    GaussLegendre1D quadrature;
    
    // Iterate over all boundary edges with this tag
    for (const BoundaryCell<1>& edge : boundaryEdges) {
        // Compute the contributions to the edge nodes using quadrature
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(edge, bc.getBoundaryFunction(t), contributions);
        
        // Add the contributions to the RHS vector
        const std::vector<unsigned int>& nodeIndices = edge.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
    }
}


template<>
void BoundaryConditions_td<1,1>::applyNeumann(
    const BoundaryCondition_td<1,1>& bc,
    const Grid<1>& mesh,
    SparseMat& /*A*/,
    VectorXd& rhs, double t)
{
    const std::vector<BoundaryCell<0>> bcs = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    if (bcs.empty()) {
        std::cerr << "Neumann 1D: no boundary cell for tag "
                  << bc.getBoundaryId() << "\n";
        return;
    }

    const BoundaryCell<0>& bc0 = bcs[0];
    const int       d = bc0.getNodeIndex(0);   // DOF globale
    const Point<1>& X = bc0.getNode(0);        // coordinata fisica

    const double g = bc.getBoundaryFunction(t).value(X); // g = μ ∂u/∂n (flusso uscente)
    rhs[d] += g;
}
