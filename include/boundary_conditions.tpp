#include "boundary_conditions.hpp"
#include <iostream>

// Costruttore con lista di condizioni
template <unsigned int dim, unsigned int returnDim>
BoundaryConditions<dim, returnDim>::BoundaryConditions(const std::vector<BoundaryCondition<dim, returnDim>>& conditions) 
    : conditions(conditions) {}


//  Adders
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::addDirichlet(int physicalTag, Function<dim, returnDim> func) {
    conditions.emplace_back(physicalTag, BCType::DIRICHLET, func);
}

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::addDirichlet(int physicalTag, Point<returnDim> value) {
    conditions.emplace_back(physicalTag, BCType::DIRICHLET, Point<returnDim>(value));
}

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::addNeumann(int physicalTag, Function<dim, returnDim> func) {
    conditions.emplace_back(physicalTag, BCType::NEUMANN, func);
}

template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::addNeumann(int physicalTag, Point<returnDim> value) {
    conditions.emplace_back(physicalTag, BCType::NEUMANN, value);
}



//Apply conditions

// Applica le condizioni al contorno al sistema lineare
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::apply(const Grid2D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
    if (conditions.empty()) {
        std::cout << "Nessuna condizione al contorno da applicare" << std::endl;
        return;
    }
    
    std::cout << "Applicazione " << conditions.size() << " condizioni al contorno..." << std::endl;

    for (const BoundaryCondition<dim, returnDim>& condition : conditions) {
        switch (condition.getType()) {
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


// Metodi helper per applicazione
template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyDirichlet(const BoundaryCondition<dim, returnDim>& bc, const Grid2D& mesh, 
                                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
    auto boundaryNodes = mesh.getBoundaryNodesByTag(bc.getPhysicalTag());
    std::cout << "  Applicando condizione di Dirichlet su tag " << bc.getPhysicalTag() 
              << " (" << boundaryNodes.size() << " nodi)" << std::endl;
    
    for (int nodeIndex : boundaryNodes) {
        // Ottieni le coordinate del nodo
        const Point<2>& nodePoint = mesh.getNode(nodeIndex);
        
        // Valuta la funzione di Dirichlet nel punto
        double dirichletValue = bc.getBoundaryFunction().value(nodePoint);

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


template <unsigned int dim, unsigned int returnDim>
void BoundaryConditions<dim, returnDim>::applyNeumann(const BoundaryCondition<dim, returnDim>& bc, const Grid2D& mesh, 
                                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {
    auto boundaryNodes = mesh.getBoundaryNodesByTag(bc.getPhysicalTag());
    std::cout << "  Applicando condizione di Neumann su tag " << bc.getPhysicalTag()
              << " (" << boundaryNodes.size() << " nodi)" << std::endl;
    
    for (int nodeIndex : boundaryNodes) {
        const Point<2>& nodePoint = mesh.getNode(nodeIndex);
        double neumannValue = bc.getBoundaryFunction().value(nodePoint);
        
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


// ---------------1D-----------------

void BoundaryConditions<1, 1>::applyDirichlet(const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
                                    Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs) {

    if (bc.getPhysicalTag() == 0){
        A.coeffRef(mesh.getStart(), mesh.getStart()) = 1.0;
        A.coeffRef(mesh.getStart(), mesh.getStart() + 1) = 0.0;
        rhs[mesh.getStart()] = bc.getBoundaryFunction().value(Point<1>(mesh.getStart()));
    } else if (bc.getPhysicalTag() == 1) {
        A.coeffRef(mesh.getEnd(), mesh.getEnd()) = 1.0;
        A.coeffRef(mesh.getEnd(), mesh.getEnd() - 1) = 0.0;
        rhs[mesh.getEnd()] = bc.getBoundaryFunction().value(Point<1>(mesh.getEnd()));
    } else {
        std::cerr << "Tag fisico non valido per condizione di Dirichlet 1D: " << bc.getPhysicalTag() << std::endl;
    }

};

void BoundaryConditions<1,1>::applyNeumann(const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs){

    if (bc.getPhysicalTag() == 0) rhs[mesh.getStart()] -= bc.getBoundaryFunction().value(Point<1>(mesh.getStart())) * 0.5;

    if (bc.getPhysicalTag() == 1) rhs[mesh.getEnd()] += bc.getBoundaryFunction().value(Point<1>(mesh.getEnd())) * 0.5;

};

// addDirichlet con Function

inline void BoundaryConditions<1,1>::addDirichlet(int physicalTag, Function<1,1> func) {
    conditions.emplace_back(physicalTag, BCType::DIRICHLET, func);
}

// addDirichlet con Point<1> (usa il costruttore specializzato di BoundaryCondition<1,1>)
inline void BoundaryConditions<1,1>::addDirichlet(int physicalTag, Point<1> value) {
    conditions.emplace_back(physicalTag, BCType::DIRICHLET, value);
}

// addNeumann con Function
inline void BoundaryConditions<1,1>::addNeumann(int physicalTag, Function<1,1> func) {
    conditions.emplace_back(physicalTag, BCType::NEUMANN, func);
}

// addNeumann con Point<1>
inline void BoundaryConditions<1,1>::addNeumann(int physicalTag, Point<1> value) {
    conditions.emplace_back(physicalTag, BCType::NEUMANN, value);
}

// apply per Grid1D (dispatch sugli helper già definiti)
inline void BoundaryConditions<1,1>::apply(const Grid1D& mesh,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs)
{
    if (conditions.empty()) {
        std::cout << "Nessuna condizione al contorno 1D\n";
        return;
    }
    for (const auto& bc : conditions) {
        if (bc.getType() == BCType::DIRICHLET)      applyDirichlet(bc, mesh, A, rhs);
        else /* NEUMANN */                          applyNeumann  (bc, mesh, A, rhs);
    }
}
