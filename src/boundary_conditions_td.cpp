#include "boundary_conditions_td.hpp"

/**
 * Apply time-dependent boundary conditions to the finite element system
 */
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions_td<dim, returnDim>::apply(const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) const {
    
    if (conditions.empty()) {
        std::cout << "No boundary condition to apply" << std::endl;
        return;
    }

    std::cout << "Applying " << conditions.size() << " boundary conditions..." << std::endl;

    // Apply Neumann first (RHS only), then Dirichlet (matrix + RHS)
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

    // Apply Dirichlet conditions last to avoid matrix overwriting
    for (const BoundaryCondition_td<dim, returnDim>& dirichletCondition : dirichletConditions) {
        applyDirichlet(dirichletCondition, mesh, A, rhs, t);
    }

    std::cout << "Boundary conditions applied successfully" << std::endl;
}

/**
 * Apply time-dependent Dirichlet boundary conditions: u = g(x,t)
 */
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions_td<dim, returnDim>::applyDirichlet(
    const BoundaryCondition_td<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) const {
    
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
        
        // Set identity: u_i = g(x_i, t)
        A.coeffRef(nodeIndex, nodeIndex) = 1.0;
        rhs[nodeIndex] = bc.getBoundaryFunction(t).value(nodePoint);
    }
}

/**
 * Apply time-dependent Neumann boundary conditions: grad(u) * n = g(x,t)
 */
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions_td<dim, returnDim>::applyNeumann(
    const BoundaryCondition_td<dim, returnDim>& bc, const Grid<dim>& mesh, 
    SparseMat& A, VectorXd& rhs, double t) const {
    
    std::vector<BoundaryCell<dim-1>> boundaryCells = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    std::cout << "  Applying Neumann boundary condition " << dim << "D on tag " << bc.getBoundaryId()
              << " (" << boundaryCells.size() << " " << (dim == 1 ? "edges" : "faces") << ")" << std::endl;

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
        for (int c = 0; c < boundaryCells.size(); ++c) {
            const BoundaryCell<dim-1>& cell = boundaryCells[c];
            
            // Integrate: ∫_∂Ω g(x,t) * φ_i dS
            std::vector<double> contributions;
            quadrature.integrateShapeFunctions(cell, bc.getBoundaryFunction(t), contributions);
            
            const std::vector<unsigned int>& nodeIndices = cell.getNodeIndexes();
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
    for (const BoundaryCell<dim-1>& cell : boundaryCells) {
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(cell, bc.getBoundaryFunction(t), contributions);
        
        const std::vector<unsigned int>& nodeIndices = cell.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
    }
#endif
}

/**
 * 1D time-dependent Neumann specialization: flux g(x,t) at boundary points
 */
template<>
void BoundaryConditions_td<1,1>::applyNeumann(
    const BoundaryCondition_td<1,1>& bc,
    const Grid<1>& mesh,
    SparseMat& /*A*/,  // Matrix A is not modified for Neumann conditions
    VectorXd& rhs, double t) const {

    const std::vector<BoundaryCell<0>> bcs = mesh.getBoundaryCellsByTag(bc.getBoundaryId());
    if (bcs.empty()) {
        std::cerr << "Neumann 1D: no boundary cell for tag "
                  << bc.getBoundaryId() << "\n";
        return;
    }

    const BoundaryCell<0>& bc0 = bcs[0];
    const int       d = bc0.getNodeIndex(0);   // Global DOF index
    const Point<1>& X = bc0.getNode(0);        // Physical coordinate

    const double g = bc.getBoundaryFunction(t).value(X);
    rhs[d] += g;
}

template class BoundaryConditions_td<1,1>;
template class BoundaryConditions_td<2,1>;
template class BoundaryConditions_td<3,1>;
