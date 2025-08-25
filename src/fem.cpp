#include "fem.hpp"

template class Fem<1>;
template class Fem<2>;
template class Fem<3>;

// Constructor
template<unsigned int dim>
Fem<dim>::Fem(Grid<dim> grid, Function<dim, 1> forcing, Function<dim, 1> diffusion, 
             Function<dim, dim> transport, Function<dim, 1> reaction,
             const BoundaryConditions<dim, 1>& boundaryConditions, QuadratureRule<dim> quadrature) :
    mesh(grid), forcing_term(forcing), diffusion_term(diffusion), 
    transport_term(transport), reaction_term(reaction),
    boundaryConditions(boundaryConditions), quadrature(quadrature)
{
#ifdef _OPENMP
    Eigen::setNbThreads(omp_get_max_threads());
#endif

    int numNodes = mesh.getNumNodes();
    A.resize(numNodes, numNodes);
    rhs.resize(numNodes);
    solution.resize(numNodes);
    rhs.setZero();
    solution.setZero();
}

// Main assembly
template<unsigned int dim>
void Fem<dim>::assemble() {
    std::cout << "### ASSEMBLE " << dim << "D ###" << std::endl;

    int numNodes = mesh.getNumNodes();
    rhs.setZero();
     
    unsigned int expNonZero = (dim==1) ? 3 : (dim==2) ? 9 : 16;
    int numElements = mesh.getNumElements();
    
    // Parallel assembly using thread-local storage
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    std::vector<std::vector<Triplet>> triplets_thread(nthreads);
    std::vector<VectorXd> rhs_thread(nthreads, VectorXd::Zero(numNodes));
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::vector<Triplet>& local_triplets = triplets_thread[tid];
        VectorXd& local_rhs = rhs_thread[tid];
        local_triplets.reserve(expNonZero * (numElements / nthreads + 1));

        #pragma omp for schedule(dynamic, 32)
        for (int e = 0; e < numElements; ++e) {
            assembleElement(e, local_triplets, local_rhs);
        }
    }
    
    // Merge RHS contributions
    for (int t = 0; t < nthreads; ++t) {
        rhs += rhs_thread[t];
    }
    
    // Merge thread-local contributions
    std::vector<Triplet> triplets;
    for (auto& v : triplets_thread) {
        triplets.insert(triplets.end(), v.begin(), v.end());
    }
#else
    std::vector<Triplet> triplets;
    triplets.reserve(expNonZero * numElements);
    for (int e = 0; e < numElements; ++e) {
        assembleElement(e, triplets, rhs);
    }
#endif
    
    A.setFromTriplets(triplets.begin(), triplets.end());
    boundaryConditions.apply(mesh, A, rhs);
}

// Element assembly
template<unsigned int dim>
void Fem<dim>::assembleElement(int elemIndex, std::vector<Triplet>& triplets, VectorXd& local_rhs) const {
    const Cell<dim>& cell = mesh.getCell(elemIndex);
    unsigned int matSize = dim+1;

    // Variables for quadrature data
    std::vector<Point<dim>> grad_phi;
    std::vector<Point<dim>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    quadrature.getQuadratureData(cell, grad_phi, quadrature_points, phi, weights);

    // Local matrices
    MatrixXd diff_local = MatrixXd::Zero(matSize, matSize); // Diffusion
    MatrixXd transport_local = MatrixXd::Zero(matSize, matSize); // Transport
    MatrixXd react_local = MatrixXd::Zero(matSize, matSize); // Reaction
    VectorXd forc_local = VectorXd::Zero(matSize); // Forcing

    // Quadrature loop
    for (int q = 0; q < quadrature_points.size(); ++q) {
        const Point<dim>& p = quadrature_points[q];
        const double w = weights[q];
        const std::vector<double>& phi_q = phi[q];
        
        // Evaluate parameters on quadrature point (cache these values)
        
        const double w_diff = diffusion_term.value(p) * w;
        const Point<dim> w_transport = transport_term.value(p) * w;
        const double w_react = reaction_term.value(p) * w;
        const double w_forc = forcing_term.value(p) * w;

        // Optimized matrix assembly with reduced memory accesses
        for (int i = 0; i < matSize; ++i) {
            const double phi_i = phi_q[i];
            const Point<dim>& grad_phi_i = grad_phi[i];
            
            forc_local(i) += w_forc * phi_i;
            
            for (int j = 0; j < matSize; ++j) {
                const double phi_j = phi_q[j];
                const Point<dim>& grad_phi_j = grad_phi[j];
                
                // Pre-compute common terms
                const double phi_i_phi_j = phi_i * phi_j;
                const double grad_dot = grad_phi_i * grad_phi_j;
                const double transport_contrib = (w_transport * grad_phi_j) * phi_i;
                
                diff_local(i,j) += w_diff * grad_dot;
                transport_local(i,j) += transport_contrib;
                react_local(i,j) += w_react * phi_i_phi_j;
            }
        }
    }
    
    // Global assembly
    for (int i = 0; i < matSize; ++i) {
        int globalI = cell.getNodeIndex(i);
        for (int j = 0; j < matSize; ++j) {
            int globalJ = cell.getNodeIndex(j);

            // Add to the global matrix only if not zero
            double value = diff_local(i,j) + transport_local(i,j) + react_local(i,j);
            if (std::abs(value) > 1e-15) {
                triplets.push_back(Triplet(globalI, globalJ, value));
            }
        }
        local_rhs[globalI] += forc_local(i);
    }
}

// Linear system solver
template<unsigned int dim>
void Fem<dim>::solve() {
    std::cout << "### SOLVE " << dim << "D ###" << std::endl;

    if(A.nonZeros() < 4e3) {
        std::cout << "[TD] Solving small system (" << A.nonZeros() << ") with SparseLU\n";
        Eigen::SparseLU<SparseMat> solver;
        solver.analyzePattern(A);
        solver.factorize(A);

        if (solver.info() != Eigen::Success) {
            std::cerr << "Factorization failed!" << std::endl;
            return;
        }
        solution = solver.solve(rhs);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solving failed!" << std::endl;
            return;
        }
    }
    else{
        // Lighter iterative solver (BiCGSTAB)
        std::cout << "[TD] Solving large system (" << A.nonZeros() << ") with BiCGSTAB\n";
        Eigen::BiCGSTAB<SparseMat, Eigen::IncompleteLUT<double>> solver;
        solver.setMaxIterations(1000);
        solver.setTolerance(1e-8);
        solver.compute(A);

        std::cout << "Finished setup for solve. Using BiCGSTAB." << std::endl;

        solution = solver.solve(rhs);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solving failed! Error is: " << solver.info() << std::endl;
            std::cerr << "Solver error at last iteration (" << solver.iterations() << "): " << solver.error() << std::endl;
            return;
        }
        else{
            std::cout << "Solving succeeded!" << std::endl;
            std::cout << "Solver iterations: " << solver.iterations() << std::endl;
            std::cout << "Solver error: " << solver.error() << std::endl;
        }
    }
}

template<unsigned int dim>
void Fem<dim>::outputCsv(const std::string& filename) const {
    if (filename.size() < 4 || filename.substr(filename.size() - 4) != ".csv") {
        std::cerr << "Error: outputCsv filename must end with '.csv'\n";
        exit(-1);
    }
    std::ofstream csvFile(filename, std::ios::out);
    if (!csvFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    csvFile << "# x " << (dim >= 2 ? "y " : "") << (dim == 3 ? "z " : "") << "u" << std::endl;
    for (unsigned int i = 0; i < mesh.getNumNodes(); ++i) {
        const Point<dim>& node = mesh.getNode(i);
        for (unsigned int j = 0; j < dim; ++j) {
            csvFile << node[j] << " ";
        }
        csvFile << solution(i) << std::endl;
    }

    std::cout << "Solution computed and written to output file" << std::endl;
    std::cout << "Solution norm: " << solution.norm() << std::endl;
}

template<unsigned int dim>
void Fem<dim>::outputVtu(const std::string& filename) const {
    if (filename.size() < 4 || filename.substr(filename.size() - 4) != ".vtu") {
        std::cerr << "Error: outputVtu filename must end with '.vtu'\n";
        exit(-1);
    }

    std::ofstream vtuFile(filename, std::ios::out);
    if (!vtuFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(-1);
    }

    vtuFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    vtuFile << "<UnstructuredGrid>\n";

    // Nodes
    vtuFile << "<Piece NumberOfPoints=\"" << mesh.getNumNodes() << "\" NumberOfCells=\"" << mesh.getNumElements() << "\">\n";
    vtuFile << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const Point<dim>& node : mesh.getUniqueNodes()) {
        for (unsigned int i = 0; i < 3; ++i) {
            vtuFile << (i<dim ? node[i] : 0.0) << " ";
        }
    }
    vtuFile << "</DataArray>\n</Points>\n";

    // Solution output as PointData
    vtuFile << "<PointData Scalars=\"solution\">\n";
    vtuFile << "<DataArray type=\"Float64\" Name=\"solution\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < mesh.getNumNodes(); ++i) {
        vtuFile << solution(i) << "\n";
    }
    vtuFile << "</DataArray>\n</PointData>\n";

    // Cells
    vtuFile << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const Cell<dim>& cell : mesh.getCells()) {
        for (unsigned int i = 0; i < cell.getN(); ++i) {
            vtuFile << cell.getNodeIndex(i) << " ";
        }
        vtuFile << "\n";
    }
    vtuFile << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    unsigned int offset = 0;
    for (const Cell<dim>& cell : mesh.getCells()) {
        offset += cell.getN();
        vtuFile << offset << "\n";
    }
    vtuFile << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    unsigned int vtkType = (dim == 1) ? 3 : // 3=line
                           (dim == 2) ? 5 : // 5=triangle
                                        10; // 10=tetrahedron
    for (unsigned int i = 0; i < mesh.getNumElements(); ++i) {
        vtuFile << vtkType << "\n";
    }
    vtuFile << "</DataArray>\n</Cells>\n";

    vtuFile << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    vtuFile.close();
}