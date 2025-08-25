#include "boundary_conditions.hpp"


// Apply all boundary conditions to the finite element system
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::apply(const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) const {
    
    if (conditions.empty()) {
        std::cout << "No boundary condition to apply" << std::endl;
        return;
    }

    std::cout << "Applying " << conditions.size() << " boundary conditions..." << std::endl;

    // Apply Neumann first (RHS only), then Dirichlet (matrix + RHS)
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

    // Apply Dirichlet conditions last to avoid matrix overwriting
    for (const BoundaryCondition<dim, returnDim>& dirichletCondition : dirichletConditions) {
        applyDirichlet(dirichletCondition, mesh, A, rhs);
    }

    std::cout << "Boundary conditions applied successfully" << std::endl;
}

// -----------------------------------------------------------------------------
// Methods for applying specific boundary condition types
// -----------------------------------------------------------------------------


// Apply Dirichlet boundary conditions: u = g on boundary

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyDirichlet(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) const {
    
    IndexVector boundaryNodes = mesh.getBoundaryNodesByTag(bc.getBoundaryId());
    std::cout << "  Applying Dirichlet condition on tag " << bc.getBoundaryId() 
              << " (" << boundaryNodes.size() << " nodes)" << std::endl;

    // Parallel safe: each thread works on different DOFs
    #pragma omp parallel for
    for (unsigned int nodeIndex : boundaryNodes) {
        const Point<dim>& nodePoint = mesh.getNode(nodeIndex);
        
        // Zero out matrix row for constrained DOF
        for (SparseMat::InnerIterator it(A, nodeIndex); it; ++it) {
            it.valueRef() = 0.0;
        }
        
        // Set identity: u_i = g(x_i)
        A.coeffRef(nodeIndex, nodeIndex) = 1.0;
        rhs[nodeIndex] = bc.getBoundaryFunction().value(nodePoint);
    }
}


// Apply Neumann boundary conditions: grad(u) * n = g on boundary

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyNeumann(
    const BoundaryCondition<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs) const {
    
    std::vector<BoundaryCell<dim-1>> boundaryFaces = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applying Neumann boundary condition " << dim << "D on tag " << bc.getBoundaryId()
              << " (" << boundaryFaces.size() << " " << (dim == 1 ? "edges" : "faces") << ")" << std::endl;

    GaussLegendre<dim-1> quadrature;
    
#ifdef _OPENMP
    // Thread-local storage to avoid race conditions
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
            
            std::vector<double> contributions;
            quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
            
            const std::vector<unsigned int>& nodeIndices = face.getNodeIndexes();
            for (size_t i = 0; i < nodeIndices.size(); ++i) {
                int globalNodeIndex = nodeIndices[i];
                rhs_local[globalNodeIndex] += contributions[i];
            }
        }
    }
    
    // Combine thread-local contributions
    for (int tid = 0; tid < nthreads; ++tid) {
        rhs += rhs_threads[tid];
    }
#else
    // Serial version: directly accumulate into global RHS
    for (const BoundaryCell<dim-1>& face : boundaryFaces) {
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(face, bc.getBoundaryFunction(), contributions);
        
        const std::vector<unsigned int>& nodeIndices = face.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
    }
#endif
}


// 1D Neumann specialization: flux condition at boundary points

template<>
void BoundaryConditions<1,1>::applyNeumann(
    const BoundaryCondition<1,1>& bc,
    const Grid<1>& mesh,
    SparseMat& /*A*/,  // Matrix A is not modified for Neumann conditions
    VectorXd& rhs) const {
        
    const std::vector<BoundaryCell<0>> bcs = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    if (bcs.empty()) {
        std::cerr << "Neumann 1D: no boundary cell for tag "
                  << bc.getBoundaryId() << "\n";
        return;
    }

    const BoundaryCell<0>& bc0 = bcs[0];
    const int       d = bc0.getNodeIndex(0);   // Global DOF index
    const Point<1>& X = bc0.getNode(0);        // Physical coordinate

    const double g = bc.getBoundaryFunction().value(X);
    rhs[d] += g;
}

template class BoundaryConditions<1,1>;
template class BoundaryConditions<2,1>;
template class BoundaryConditions<3,1>;
