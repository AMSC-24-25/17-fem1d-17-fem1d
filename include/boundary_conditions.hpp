#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "function.hpp"
#include "point.hpp"
#include "grid2D.hpp"
#include <Eigen/Sparse>
#include <vector>

// Enum per il tipo di condizione al contorno
enum class BCType {
    DIRICHLET,  // u = valore specificato
    NEUMANN     // du/dn = valore specificato
};

// Configurazione per una singola condizione al contorno
struct BoundaryCondition {
    int physicalTag;
    BCType type;
    Function<2> boundaryFunction;
    
    // Costruttore per funzione
    BoundaryCondition(int tag, BCType bcType, Function<2> func) 
        : physicalTag(tag), type(bcType), boundaryFunction(func) {}
    
    // Costruttore per valore costante (per comodità)
    BoundaryCondition(int tag, BCType bcType, double value) 
        : physicalTag(tag), type(bcType), boundaryFunction([value](Point<2> p) { return value; }) {}
};

// Classe per gestire tutte le condizioni al contorno
class BoundaryConditions {
public:
    // Costruttori
    BoundaryConditions() = default;
    BoundaryConditions(const std::vector<BoundaryCondition>& conditions);
    
    // Applica le condizioni al contorno al sistema lineare
    void apply(const Grid2D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    
    // Aggiungi condizioni
    void addDirichlet(int physicalTag, Function<2> func);
    void addDirichlet(int physicalTag, double value);
    void addNeumann(int physicalTag, Function<2> func);
    void addNeumann(int physicalTag, double value);

private:
    std::vector<BoundaryCondition> conditions;
    
    // Metodi helper per applicazione
    void applyDirichlet(const BoundaryCondition& bc, const Grid2D& mesh, 
                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    void applyNeumann(const BoundaryCondition& bc, const Grid2D& mesh, 
                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
};

#endif // BOUNDARY_CONDITIONS_HPP
