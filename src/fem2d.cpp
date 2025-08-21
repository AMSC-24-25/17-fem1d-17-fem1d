#include "fem2d.hpp"

template class Fem<2>;

// Costruttore moderno con BoundaryConditions
template<unsigned int dim>
Fem<dim>::Fem(Grid<dim> grid, Function<dim, 1> forcing, Function<dim, 1> diffusion, 
             Function<dim, dim> transport, Function<dim, 1> reaction,
             const BoundaryConditions<dim, 1>& boundaryConditions) :
    mesh(grid), forcing_term(forcing), diffusion_term(diffusion), 
    transport_term(transport), reaction_term(reaction),
    boundaryConditions(boundaryConditions)
{
    // Inizializza matrici
    int numNodes = mesh.getNumNodes();
    A.resize(numNodes, numNodes);
    rhs.resize(numNodes);
    solution.resize(numNodes);
    rhs.setZero();
    solution.setZero();
}

template<>
std::unique_ptr<QuadratureRule<2>> Fem<2>::getQuadratureRule() const {
    return std::make_unique<OrderTwoQuadrature<2>>();
}

// Main assembly
template<unsigned int dim>
void Fem<dim>::assemble() {
    std::cout << "### ASSEMBLE " << dim << "D ###" << std::endl;

    // Inizializzazione
    int numNodes = mesh.getNumNodes();
    rhs.setZero();
    
    std::vector<Triplet> triplets;
    unsigned int expNonZero = (dim==1) ? 3 : (dim==2) ? 9 : 25; // Estimate non-zero entries per element
    triplets.reserve(expNonZero * mesh.getNumElements());

    // Obtain quadrature rule from specialised method
    std::unique_ptr<QuadratureRule<dim>> quadrature = getQuadratureRule();

    // Loop on elements - assemble local matrix for each. Appends non-zero values to triplets
    for (unsigned int e = 0; e < mesh.getNumElements(); ++e) {
        assembleElement(e, *quadrature, triplets);
    }
    
    // Global sparse matrix assembly
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Apply boundary conditions
    boundaryConditions.apply(mesh, A, rhs);
    
    std::cout << "Matrix assembled: " << A.rows() << "x" << A.cols() 
              << " with " << A.nonZeros() << " non-zeros" << std::endl;
}

// Assemblaggio di un singolo elemento (triangolo)
template<unsigned int dim>
void Fem<dim>::assembleElement(int elemIndex, QuadratureRule<dim>& quad, 
                           std::vector<Triplet>& triplets) {

    const Cell<dim>& cell = mesh.getCell(elemIndex);

    // // Verifica che sia effettivamente un triangolo
    // if (cell.getN() != 3) {
    //     std::cerr << "ERROR: Element " << elemIndex << " is not a triangle (has " 
    //               << cell.getN() << " nodes)" << std::endl;
    //     return; // Skip this element
    // }
    
    unsigned int matSize = (dim == 2) ? 3 : 5;

    // Variables for quadrature data
    std::vector<Point<dim>> grad_phi;
    std::vector<Point<dim>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    // Get quadrature data
    quad.getQuadratureData(cell, grad_phi, quadrature_points, phi, weights);

    // Local matrices for element
    MatrixXd diff_local = MatrixXd::Zero(matSize, matSize); // Diffusione
    MatrixXd transport_local = MatrixXd::Zero(matSize, matSize); // Massa/Reazione
    MatrixXd react_local = MatrixXd::Zero(matSize, matSize); // Massa/Reazione
    VectorXd forc_local = VectorXd::Zero(matSize); // Termine forzante

    // Loop on quadrature points
    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        const Point<dim>& p = quadrature_points[q];
        double w = weights[q];
        
        // Evaluate parameters on quadrature point
        double diff_val = diffusion_term.value(p);
        double react_val = reaction_term.value(p);
        double forc_val = forcing_term.value(p);

        // Contributions to local matrices
        for (int i = 0; i < matSize; ++i) {
            for (int j = 0; j < matSize; ++j) {
                
                // Diffusion contribution matrix: ∫ ∇φᵢ · ∇φⱼ dx
                // Transport contribution term ∫ (b·∇φᵢ) φⱼ dx
                // Note: product of 2 Points is their scalar product
                diff_local(i,j) += w * diff_val * (grad_phi[i] * grad_phi[j]);
                transport_local(i,j) += w * (transport_term.value(p) * grad_phi[i]) * phi[q][j];

                // Reaction contribution matrix: ∫ φᵢ φⱼ dx
                react_local(i,j) += w * react_val * phi[q][i] * phi[q][j];
            }
            // Forcing term: ∫ f φᵢ dx
            forc_local(i) += w * forc_val * phi[q][i];
        }
    }
    
    // Assembly in global matrix
    for (int i = 0; i < matSize; ++i) {
        for (int j = 0; j < matSize; ++j) {
            int globalI = cell.getNodeIndex(i);
            int globalJ = cell.getNodeIndex(j);

            // Aggiungi alla matrice globale solo se non zero
            double value = diff_local(i,j) + react_local(i,j);
            if (std::abs(value) > 1e-14) {
                triplets.push_back(Triplet(globalI, globalJ, value));
            }
        }

        // Assembly of RHS
        int globalI = cell.getNodeIndex(i);
        rhs[globalI] += forc_local(i);
    }
}

// Risoluzione del sistema lineare
template<unsigned int dim>
void Fem<dim>::solve() {
    std::cout << "### SOLVE " << dim << "D ###" << std::endl;

    // Risolvi il sistema Ax = b usando SparseLU
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