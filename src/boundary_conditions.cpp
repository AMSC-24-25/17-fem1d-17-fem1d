#include "boundary_conditions.hpp"

// -----------------------------------------------------------------------------
// Apply boundary conditions (main method)
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::apply(const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) const {
    
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
// Methods for applying specific boundary condition types
// -----------------------------------------------------------------------------

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyDirichlet(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) const {
    
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

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyNeumann(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) const {
    
    // Get boundary faces with the specified physical tag
    std::vector<BoundaryCell<dim-1>> boundaryFaces = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applying Neumann boundary condition " << dim << "D on tag " << bc.getBoundaryId()
              << " (" << boundaryFaces.size() << " " << (dim == 1 ? "edges" : "faces") << ")" << std::endl;

    // Initialize 2D quadrature for triangular faces
    GaussLegendre<dim-1> quadrature;
    
    // Iterate over all boundary faces with this tag
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    static std::vector<VectorXd> rhs_threads(nthreads, VectorXd::Zero(rhs.size()));
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        VectorXd& rhs_local = rhs_threads[tid];
        rhs_local.setZero();
        
        #pragma omp for schedule(static)
        for (int f = 0; f < boundaryFaces.size(); ++f) {
            const BoundaryCell<dim-1>& face = boundaryFaces[f];
            
            // Compute the contributions to the face nodes using quadrature
            std::vector<double> contributions;
            quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
            
            // Add the contributions to the local RHS vector
            const std::vector<unsigned int>& nodeIndices = face.getNodeIndexes();
            for (size_t i = 0; i < nodeIndices.size(); ++i) {
                int globalNodeIndex = nodeIndices[i];
                rhs_local[globalNodeIndex] += contributions[i];
            }
        }
    }
    
    // Sum all thread-local vectors
    for (int tid = 0; tid < nthreads; ++tid) {
        rhs += rhs_threads[tid];
    }
#else
    for (const BoundaryCell<dim-1>& face : boundaryFaces) {
        // Compute the contributions to the face nodes using quadrature
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
        
        // Add the contributions to the RHS vector
        const std::vector<unsigned int>& nodeIndices = face.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
    }
#endif
}

template<>
void BoundaryConditions<1,1>::applyNeumann(
    const BoundaryCondition<1,1>& bc,
    const Grid<1>& mesh,
    SparseMat& /*A*/,
    VectorXd& rhs) const {
        
    const std::vector<BoundaryCell<0>> bcs = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    if (bcs.empty()) {
        std::cerr << "Neumann 1D: nessun boundary cell per tag "
                  << bc.getBoundaryId() << "\n";
        return;
    }

    const BoundaryCell<0>& bc0 = bcs[0];
    const int       d = bc0.getNodeIndex(0);   // DOF globale
    const Point<1>& X = bc0.getNode(0);        // coordinata fisica

    const double g = bc.getBoundaryFunction().value(X); // g = μ ∂u/∂n (flusso uscente)
    rhs[d] += g;
}

template class BoundaryConditions<1,1>;
template class BoundaryConditions<2,1>;
template class BoundaryConditions<3,1>;
