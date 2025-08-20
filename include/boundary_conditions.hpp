#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "function.hpp"
#include "point.hpp"
#include "grid2D.hpp"
#include "grid1D.hpp"
#include <Eigen/Sparse>
#include <vector>

// Enum per il tipo di condizione al contorno
enum class BCType {
    DIRICHLET,  // u = valore specificato
    NEUMANN     // du/dn = valore specificato
};


template<unsigned int dim, unsigned int returnDim>
class BoundaryCondition {

private:
    int physicalTag;
    BCType type;
    Function<dim, returnDim> boundaryFunction;

public:
    // Costruttore per funzione
    BoundaryCondition(int tag, BCType bcType, Function<dim, returnDim> func) 
        : physicalTag(tag), type(bcType), boundaryFunction(func) {}
    
    // Costruttore per valore costante (per comodità)
    BoundaryCondition(int tag, BCType bcType, Point<returnDim> value) 
        : physicalTag(tag), type(bcType), boundaryFunction([value](Point<dim> p) { return value; }) {}

    int getPhysicalTag() const { return physicalTag; }
    BCType getType() const { return type; }
    Function<dim, returnDim> getBoundaryFunction() const { return boundaryFunction; }
};



template <unsigned int dim, unsigned int returnDim>
class BoundaryConditions {
public:
    // Costruttori
    BoundaryConditions() = default;
    BoundaryConditions(const std::vector<BoundaryCondition<dim, returnDim>>& conditions);

    // Applica le condizioni al contorno al sistema lineare
    void apply(const Grid2D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    
    // Aggiungi condizioni
    void addDirichlet(int physicalTag, Function<dim, returnDim> func);
    void addDirichlet(int physicalTag, Point<returnDim> value);
    void addNeumann(int physicalTag, Function<dim, returnDim> func);
    void addNeumann(int physicalTag, Point<returnDim> value);

private:
    std::vector<BoundaryCondition<dim, returnDim>> conditions;

    // Metodi helper per applicazione
    void applyDirichlet(const BoundaryCondition<dim, returnDim>& bc, const Grid2D& mesh, 
                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    void applyNeumann(const BoundaryCondition<dim, returnDim>& bc, const Grid2D& mesh, 
                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
};


template<>
class BoundaryConditions <1, 1> {
public:
    // Costruttori
    BoundaryConditions() = default;
    BoundaryConditions(const std::vector<BoundaryCondition<1, 1>>& conditions) 
        : conditions(conditions) {}
    // Applica le condizioni al contorno al sistema lineare
    void apply(const Grid1D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    
    // Aggiungi condizioni
    void addDirichlet(int physicalTag, Function<1, 1> func);
    void addDirichlet(int physicalTag, Point<1> value);
    void addNeumann(int physicalTag, Function<1, 1> func);
    void addNeumann(int physicalTag, Point<1> value);

private:
    std::vector<BoundaryCondition<1, 1>> conditions;

    // Metodi helper per applicazione
    inline void applyDirichlet(const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    inline void applyNeumann(const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
};


#include "boundary_conditions.tpp"

#endif // BOUNDARY_CONDITIONS_HPP
