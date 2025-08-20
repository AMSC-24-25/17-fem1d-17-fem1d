#include "../include/fem2d.hpp"

// Costruttore moderno con BoundaryConditions
Fem2D::Fem2D(Grid2D grid, Function<2> forcing, Function<2> diffusion, 
             Function<2> transport, Function<2> reaction,
             const BoundaryConditions& boundaryConditions) :
    mesh(grid), forcing_term(forcing), diffusion_term(diffusion), 
    transport_term(transport), reaction_term(reaction),
    boundaryConditions(boundaryConditions),
    isNeumann1(false), isNeumann2(false), 
    boundary1([](Point<2> p) { return 0.0; }), 
    boundary2([](Point<2> p) { return 0.0; })
{
    // Inizializza matrici
    int numNodes = mesh.getNumNodes();
    A.resize(numNodes, numNodes);
    rhs.resize(numNodes);
    solution.resize(numNodes);
    rhs.setZero();
    solution.setZero();
}

// Costruttore legacy (per compatibilità)
Fem2D::Fem2D(Grid2D grid, Function<2> forcing, Function<2> diffusion, 
             Function<2> transport, Function<2> reaction,
             bool isNeumann1, bool isNeumann2, 
             Function<2> boundary1, Function<2> boundary2) :
    mesh(grid), forcing_term(forcing), diffusion_term(diffusion), 
    transport_term(transport), reaction_term(reaction),
    isNeumann1(isNeumann1), isNeumann2(isNeumann2), 
    boundary1(boundary1), boundary2(boundary2)
{
    // Inizializza matrici
    int numNodes = mesh.getNumNodes();
    A.resize(numNodes, numNodes);
    rhs.resize(numNodes);
    solution.resize(numNodes);
    rhs.setZero();
    solution.setZero();
}

// Assemblaggio principale
void Fem2D::assemble() {
    std::cout << "### ASSEMBLE 2D ###" << std::endl;
    
    // Inizializzazione
    int numNodes = mesh.getNumNodes();
    rhs.setZero();
    
    std::vector<Triplet> triplets;
    triplets.reserve(9 * mesh.getNumElements()); // Stima: 9 entries per triangolo
    
    // Regola di quadratura per triangoli
    orderTwoQuadrature quadrature;
    
    // Loop sui triangoli - CUORE DELL'ASSEMBLAGGIO
    for (unsigned int e = 0; e < mesh.getNumElements(); ++e) {
        assembleElement(e, quadrature, triplets);
    }
    
    // Assemblaggio finale della matrice sparsa
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    // Applica condizioni al contorno
    boundaryConditions.apply(mesh, A, rhs);
    
    std::cout << "Matrix assembled: " << A.rows() << "x" << A.cols() 
              << " with " << A.nonZeros() << " non-zeros" << std::endl;
}

// Assemblaggio di un singolo elemento (triangolo)
void Fem2D::assembleElement(int elemIndex, BarycentricQuadRule& quad, 
                           std::vector<Triplet>& triplets) {
    
    const Cell<2>& triangle = mesh.getCell(elemIndex);
    
    // Verifica che sia effettivamente un triangolo
    if (triangle.getN() != 3) {
        std::cerr << "ERROR: Element " << elemIndex << " is not a triangle (has " 
                  << triangle.getN() << " nodes)" << std::endl;
        return; // Skip this element
    }
    
    // Calcola matrici locali 3x3 usando la quadratura esistente
    
    // Variabili per i dati di quadratura
    std::vector<Point<2>> grad_phi;
    std::vector<Point<2>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    // Ottieni dati di quadratura
    quad.getQuadratureData(triangle, grad_phi, quadrature_points, phi, weights);
    
    // Matrici locali 3x3 per triangolo P1
    Matrix3d K_local = Matrix3d::Zero(); // Diffusione
    Matrix3d M_local = Matrix3d::Zero(); // Massa/Reazione
    Vector3d f_local = Vector3d::Zero(); // Termine forzante
    
    // Loop sui punti di quadratura
    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        const Point<2>& p = quadrature_points[q];
        double w = weights[q];
        
        // Valutazioni delle funzioni nel punto di quadratura
        double diff_val = diffusion_term.value(p);
        double react_val = reaction_term.value(p);
        double forc_val = forcing_term.value(p);
        
        // Contributi alle matrici locali
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // Matrice di rigidezza: ∫ ∇φᵢ · ∇φⱼ dx
                K_local(i,j) += w * diff_val * (grad_phi[i][0] * grad_phi[j][0] + 
                                              grad_phi[i][1] * grad_phi[j][1]);
                // Matrice di massa: ∫ φᵢ φⱼ dx  
                M_local(i,j) += w * react_val * phi[q][i] * phi[q][j];
            }
            // Termine forzante: ∫ f φᵢ dx
            f_local(i) += w * forc_val * phi[q][i];
        }
    }
    
    // Assemblaggio nella matrice globale
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int globalI = triangle.getNodeIndex(i);
            int globalJ = triangle.getNodeIndex(j);

            // Aggiungi alla matrice globale solo se non zero
            double value = K_local(i,j) + M_local(i,j);
            if (std::abs(value) > 1e-14) {
                triplets.push_back(Triplet(globalI, globalJ, value));
            }
        }
        
        // Assemblaggio del termine noto (RHS)
        int globalI = triangle.getNodeIndex(i);
        rhs[globalI] += f_local(i);
    }
}

// Risoluzione del sistema lineare
void Fem2D::solve(std::ofstream& output) {
    std::cout << "### SOLVE 2D ###" << std::endl;
    
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
    
    // Scrivi soluzione su file (formato per visualizzazione)
    output << "# x y u" << std::endl;
    for (unsigned int i = 0; i < mesh.getNumNodes(); ++i) {
        const Point<2>& node = mesh.getNode(i);
        output << node[0] << " " << node[1] << " " << solution(i) << std::endl;
    }
    
    std::cout << "Solution computed and written to output file" << std::endl;
    std::cout << "Solution norm: " << solution.norm() << std::endl;
}

void Fem2D::outputVtk(const std::string& filename) const {
    if (filename.size() < 4 || filename.substr(filename.size() - 4) != ".vtu") {
        std::cerr << "Error: outputVtk filename must end with '.vtu'\n";
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
    for (const Point<2>& node : mesh.getUniqueNodes()) {
        vtuFile << node.x() << " " << node.y() << " 0.0\n"; // Z coordinate required 
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
    for (const Cell<2>& cell : mesh.getCells()) {
        for (unsigned int i = 0; i < cell.getN(); ++i) {
            vtuFile << cell.getNodeIndex(i) << " ";
        }
        vtuFile << "\n";
    }
    vtuFile << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    unsigned int offset = 0;
    for (const Cell<2>& cell : mesh.getCells()) {
        offset += cell.getN();
        vtuFile << offset << "\n";
    }
    vtuFile << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < mesh.getNumElements(); ++i) {
        vtuFile << "5\n"; // VTK_TRIANGLE
    }
    vtuFile << "</DataArray>\n</Cells>\n";

    vtuFile << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    vtuFile.close();
}