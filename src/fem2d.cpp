#include "../include/fem2d.hpp"

// Costruttore
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
    
    std::cout << "Fem2D initialized with " << numNodes << " nodes and " 
              << mesh.getNumElements() << " elements" << std::endl;
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
    BarycentricQuadRule quad = quadTriOrder2();
    
    // Loop sui triangoli - CUORE DELL'ASSEMBLAGGIO
    for (unsigned int e = 0; e < mesh.getNumElements(); ++e) {
        assembleElement(e, quad, triplets);
    }
    
    // Assemblaggio finale della matrice sparsa
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    // Applica condizioni al contorno
    applyBoundaryConditions();
    
    std::cout << "Matrix assembled: " << A.rows() << "x" << A.cols() 
              << " with " << A.nonZeros() << " non-zeros" << std::endl;
}

// Assemblaggio di un singolo elemento (triangolo)
void Fem2D::assembleElement(int elemIndex, BarycentricQuadRule& quad, 
                           std::vector<Triplet>& triplets) {
    
    const Cell<2>& triangle = mesh.getCell(elemIndex);
    
    // Calcola matrici locali 3x3 usando la quadratura esistente
    Matrix3d diffusionLocal, transportLocal, reactionLocal;
    Vector3d forcingLocal;
    
    // Converte transport_term in std::function per compatibilità
    auto transportFunc = [this](const Point<2>& p) -> Point<2> {
        double val = transport_term(p);
        return Point<2>(val, 0.0); // Assume trasporto solo in direzione x
    };
    
    // Calcola le matrici locali usando la quadratura
    quad.localMatricesP1(
        triangle,
        diffusion_term,    // Function<2>
        reaction_term,     // Function<2>  
        transportFunc,     // std::function<Point<2>(Point<2>)>
        forcing_term,      // Function<2>
        diffusionLocal,    // Output: matrice 3x3
        transportLocal,    // Output: matrice 3x3  
        reactionLocal,     // Output: matrice 3x3
        forcingLocal       // Output: vettore 3x1
    );
    
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

// Applica condizioni al contorno
void Fem2D::applyBoundaryConditions() {
    std::cerr << "TODO Boundary Conditions" << std::endl;
    
    // TODO: Implementare condizioni al contorno 2D
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
