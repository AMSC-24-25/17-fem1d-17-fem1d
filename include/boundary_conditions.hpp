#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include "function.hpp"
#include "point.hpp"
#include "grid1D.hpp"
#include "grid.hpp"
#include <Eigen/Sparse>
#include <vector>

// =============================================================================
// ENUMS AND BASIC STRUCTURES
// =============================================================================

/// Enum per il tipo di condizione al contorno
enum class BCType {
    DIRICHLET,  ///< u = valore specificato
    NEUMANN     ///< du/dn = valore specificato
};

// =============================================================================
// BOUNDARY CONDITION CLASS (GENERIC)
// =============================================================================

/// Rappresenta una singola condizione al contorno
template<unsigned int dim, unsigned int returnDim>
class BoundaryCondition {
private:
    int physicalTag;                            ///< Tag fisico del bordo
    BCType type;                               ///< Tipo di condizione (Dirichlet/Neumann)
    Function<dim, returnDim> boundaryFunction; ///< Funzione che definisce la condizione

public:
    /// Costruttore per funzione
    BoundaryCondition(int tag, BCType bcType, Function<dim, returnDim> func) 
        : physicalTag(tag), type(bcType), boundaryFunction(func) {}
    
    /// Costruttore per valore costante (per comodità)
    BoundaryCondition(int tag, BCType bcType, Point<returnDim> value) 
        : physicalTag(tag), type(bcType), 
          boundaryFunction([value](Point<dim> p) { return value; }) {}

    // Getters
    int getPhysicalTag() const { return physicalTag; }
    BCType getType() const { return type; }
    Function<dim, returnDim> getBoundaryFunction() const { return boundaryFunction; }
};

// =============================================================================
// BOUNDARY CONDITIONS CLASS (GENERIC - FOR 2D/3D)
// =============================================================================

/// Gestisce un insieme di condizioni al contorno per problemi 2D/3D
template <unsigned int dim, unsigned int returnDim>
class BoundaryConditions {
public:
    // Costruttori
    BoundaryConditions() = default;
    BoundaryConditions(const std::vector<BoundaryCondition<dim, returnDim>>& conditions);

    // Metodi per aggiungere condizioni
    void addDirichlet(int physicalTag, Function<dim, returnDim> func);
    void addDirichlet(int physicalTag, Point<returnDim> value);
    void addNeumann(int physicalTag, Function<dim, returnDim> func);
    void addNeumann(int physicalTag, Point<returnDim> value);

    // Applicazione delle condizioni al contorno
    void apply(const Grid2D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);

private:
    std::vector<BoundaryCondition<dim, returnDim>> conditions;

    // Metodi helper per applicazione
    void applyDirichlet(const BoundaryCondition<dim, returnDim>& bc, const Grid2D& mesh, 
                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    void applyNeumann(const BoundaryCondition<dim, returnDim>& bc, const Grid2D& mesh, 
                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
};

// =============================================================================
// BOUNDARY CONDITIONS CLASS (SPECIALIZATION FOR 1D)
// =============================================================================

/// Specializzazione per problemi 1D
template<>
class BoundaryConditions<1, 1> {
public:
    // Costruttori
    BoundaryConditions() = default;
    BoundaryConditions(const std::vector<BoundaryCondition<1, 1>>& conditions) 
        : conditions(conditions) {}

    // Metodi per aggiungere condizioni
    inline void addDirichlet(int physicalTag, Function<1, 1> func);
    inline void addDirichlet(int physicalTag, Point<1> value);
    inline void addNeumann(int physicalTag, Function<1, 1> func);
    inline void addNeumann(int physicalTag, Point<1> value);

    // Applicazione delle condizioni al contorno
    inline void apply(const Grid1D& mesh, Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);

private:
    std::vector<BoundaryCondition<1, 1>> conditions;

    // Metodi helper per applicazione
    inline void applyDirichlet(const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
                       Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
    inline void applyNeumann(const BoundaryCondition<1, 1>& bc, const Grid1D& mesh, 
                     Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& rhs);
};

// =============================================================================
// TEMPLATE IMPLEMENTATION
// =============================================================================

#include "boundary_conditions.tpp"

#endif // BOUNDARY_CONDITIONS_HPP
