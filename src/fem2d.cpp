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
    orderTwoQuadrature quadrature = orderTwoQuadrature();

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
void Fem2D::assembleElement(int elemIndex, orderTwoQuadrature& quadrature, 
                           std::vector<Triplet>& triplets) {
    
    const Cell<2>& triangle = mesh.getCell(elemIndex);
    
    // Verifica che sia effettivamente un triangolo
    if (triangle.getN() != 3) {
        std::cerr << "ERROR: Element " << elemIndex << " is not a triangle (has " 
                  << triangle.getN() << " nodes)" << std::endl;
        return; // Skip this element
    }
    
    // Calcola matrici locali 3x3 usando la quadratura esistente
    Matrix3d diffusionLocal, transportLocal, reactionLocal;
    Vector3d forcingLocal;

    diffusionLocal.setZero(); transportLocal.setZero(); 
    reactionLocal.setZero(); forcingLocal.setZero();
    
    //TODO
    // Converte transport_term in std::function per compatibilità
    auto transportFunc = [this](const Point<2>& p) -> Point<2> {
        double val = transport_term(p);
        return Point<2>(val, 0.0); // Assume trasporto solo in direzione x
    };

    std::vector<Point<2>> grad_phi;
    std::vector<Point<2>> quadrature_points;
    std::vector<std::vector<double>> phi_vector;
    std::vector<double> weights;

    // Calcola le matrici locali usando la quadratura
    quadrature.getQuadratureData(
        triangle,
        grad_phi,
        quadrature_points,
        phi_vector,
        weights
    );
    for (int q = 0; q < quadrature_points.size(); ++q) {
        const Point<2>& p = quadrature_points[q];
        double weight = weights[q];
        double diff_loc = diffusion_term(p);
        double react_loc = reaction_term(p);
        Point<2> b = transportFunc(p);
        std::vector<double> phi = phi_vector[q];
        for(int i=0;i<3;++i){
            double phi_i = phi[i];
            //LHS
            for(int j=0;j<3;++j){
                double phi_j = phi[j];
                double gradDot = (grad_phi[j][0]*grad_phi[i][0] + grad_phi[j][1]*grad_phi[i][1]);
                double bDot = (b[0]*grad_phi[j][0] + b[1]*grad_phi[j][1]);
                diffusionLocal(i,j) += diff_loc * gradDot * weight;
                transportLocal(i,j) += bDot * phi_i * weight; // check sign via weak form
                reactionLocal(i,j) += react_loc * phi_i * phi_j * weight;
            }
            
            // RHS
            forcingLocal[i] += forcing_term(p)*phi_i*weight;
        }
    }

    // Matrice locale totale: A_local = K + M + C
    Matrix3d A_local = diffusionLocal + transportLocal + reactionLocal;
    
    // Assemblaggio nella matrice globale
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int globalI = mesh.getCell(elemIndex).getNodeIndex(i);
            int globalJ = mesh.getCell(elemIndex).getNodeIndex(j);

            // Aggiungi alla matrice globale solo se non zero
            if (std::abs(A_local(i,j)) > 1e-14) {
                triplets.push_back(Triplet(globalI, globalJ, A_local(i,j)));
            }
        }
        
        // Assemblaggio del termine noto (RHS)
        int globalI = mesh.getCell(elemIndex).getNodeIndex(i);
        rhs[globalI] += forcingLocal[i];
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
