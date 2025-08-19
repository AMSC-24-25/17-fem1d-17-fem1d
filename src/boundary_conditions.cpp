#include "../include/boundary_conditions.hpp"
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
    auto boundaryNodes = mesh.getBoundaryNodesByTag(bc.physicalTag);
    std::cout << "  Applicando condizione di Neumann su tag " << bc.physicalTag 
              << " (" << boundaryNodes.size() << " nodi)" << std::endl;
    
    for (int nodeIndex : boundaryNodes) {
        const Point<2>& nodePoint = mesh.getNode(nodeIndex);
        double neumannValue = bc.boundaryFunction.value(nodePoint);
        
        // TODO: Implementare correttamente le condizioni di Neumann
        // Richiedono integrazione sui boundary elements e calcolo delle derivate normali
        std::cout << "    Neumann BC at node " << nodeIndex 
                 << " with flux " << neumannValue << " (placeholder - non implementato)" << std::endl;
        
        // Per ora non modifichiamo la matrice A per le condizioni di Neumann
        // L'implementazione completa richiederebbe:
        // 1. Identificazione degli elementi di bordo
        // 2. Calcolo delle derivate normali
        // 3. Integrazione sui bordi
        // 4. Aggiunta del contributo al RHS
    }
}
