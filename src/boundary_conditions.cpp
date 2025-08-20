#include "../include/boundary_conditions.hpp"
#include "../include/quadrature.hpp"
#include <iostream>

// Costruttore con lista di condizioni
BoundaryConditions::BoundaryConditions(const std::vector<BoundaryCondition>& conditions) 
    : conditions(conditions) {
}

// Applica le condizioni al contorno al sistema lineare
void BoundaryConditions::apply(const Grid2D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
    if (conditions.empty()) {
        std::cout << "Nessuna condizione al contorno da applicare" << std::endl;
        return;
    }
    
    std::cout << "Applicazione " << conditions.size() << " condizioni al contorno..." << std::endl;
    
    for (const auto& condition : conditions) {
        switch (condition.type) {
            case BCType::DIRICHLET:
                applyDirichlet(condition, mesh, A, rhs);
                break;
            case BCType::NEUMANN:
                applyNeumann(condition, mesh, A, rhs);
                break;
        }
    }
    
    std::cout << "Condizioni al contorno applicate con successo" << std::endl;
}

// Aggiungi condizioni
void BoundaryConditions::addDirichlet(int physicalTag, Function<2> func) {
    conditions.emplace_back(physicalTag, BCType::DIRICHLET, func);
}

void BoundaryConditions::addDirichlet(int physicalTag, double value) {
    conditions.emplace_back(physicalTag, BCType::DIRICHLET, value);
}

void BoundaryConditions::addNeumann(int physicalTag, Function<2> func) {
    conditions.emplace_back(physicalTag, BCType::NEUMANN, func);
}

void BoundaryConditions::addNeumann(int physicalTag, double value) {
    conditions.emplace_back(physicalTag, BCType::NEUMANN, value);
}

// Metodi helper per applicazione
void BoundaryConditions::applyDirichlet(const BoundaryCondition& bc, const Grid2D& mesh, 
                                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
    auto boundaryNodes = mesh.getBoundaryNodesByTag(bc.physicalTag);
    std::cout << "  Applicando condizione di Dirichlet su tag " << bc.physicalTag 
              << " (" << boundaryNodes.size() << " nodi)" << std::endl;
    
    for (int nodeIndex : boundaryNodes) {
        // Ottieni le coordinate del nodo
        const Point<2>& nodePoint = mesh.getNode(nodeIndex);
        
        // Valuta la funzione di Dirichlet nel punto
        double dirichletValue = bc.boundaryFunction.value(nodePoint);
        
        // Applica condizione di Dirichlet: u[nodeIndex] = dirichletValue
        
        // Azzeramento della riga
        for (int k = 0; k < A.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, k); it; ++it) {
                if (it.row() == nodeIndex) {
                    it.valueRef() = 0.0;
                }
            }
        }
        
        // Imposta 1 sulla diagonale
        A.coeffRef(nodeIndex, nodeIndex) = 1.0;
        
        // Imposta il valore calcolato nel RHS
        rhs[nodeIndex] = dirichletValue;
    }
}

void BoundaryConditions::applyNeumann(const BoundaryCondition& bc, const Grid2D& mesh, 
                                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
    // Ottieni gli edge di bordo con il tag fisico specifico
    auto boundaryEdges = mesh.getBoundaryEdgesByTag(bc.physicalTag);
    std::cout << "  Applicando condizione di Neumann su tag " << bc.physicalTag 
              << " (" << boundaryEdges.size() << " edge)" << std::endl;
    
    // Inizializza la quadratura 1D
    GaussLegendre1D quadrature;
    
    // Itera su tutti gli edge di bordo con questo tag
    for (const auto& edge : boundaryEdges) {
        // Calcola i contributi ai nodi dell'edge usando la quadratura
        std::vector<double> contributions;
        quadrature.integrateShapeFunctions(edge, bc.boundaryFunction, contributions);
        
        // Aggiungi i contributi al vettore RHS
        const auto& nodeIndices = edge.getNodeIndexes();
        for (size_t i = 0; i < nodeIndices.size(); ++i) {
            int globalNodeIndex = nodeIndices[i];
            rhs[globalNodeIndex] += contributions[i];
        }
        
        std::cout << "    Edge con nodi [" << nodeIndices[0] << ", " << nodeIndices[1] 
                 << "] - contributi: [" << contributions[0] << ", " << contributions[1] << "]" << std::endl;
    }
}

// TODO: Implementare metodi per 1D e 3D
// Per 1D: applica condizioni di Neumann direttamente ai punti di bordo (già implementato in Fem1D)
// void applyNeumann1D(const BoundaryCondition& bc, const Grid1D& mesh, 
//                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);

// TODO: Per 3D: applica condizioni di Neumann integrando sulle facce di bordo
// void applyNeumann3D(const BoundaryCondition& bc, const Grid3D& mesh, 
//                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
//     // 1. Ottenere le facce di bordo con getBoundaryFacesByTag(bc.physicalTag)
//     // 2. Usare quadratura 2D per integrare ∫ g(x)·φᵢ(x) dS su ogni faccia
//     // 3. Assemblare i contributi nel RHS come fatto per il 2D
//     // Richiede: GaussLegendre2D quadrature, BoundaryCell<2> per facce, 
//     //          mapToGlobalFace(), getShapeFunctions2D(), faceArea()
// }
