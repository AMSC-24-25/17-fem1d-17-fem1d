#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>
#include "grid1D.hpp"
#include "grid.hpp"
#include "function.hpp"
#include "boundary_conditions.hpp"
#include <memory>

// Struttura per la configurazione del problema
struct ProblemConfig {
    unsigned int dimension;
    std::string mesh_file;
    std::string output_file;
    int grid_size;  // For 1D uniform grids
};

// Struttura per la configurazione dell'equazione
struct EquationConfig {
    // Unified approach: everything is a function expression
    std::string diffusion_function;
    std::string transport_function;
    std::string reaction_function;
    std::string forcing_function;
    
    // Scalar coefficients for backward compatibility
    double diffusion_coefficient = 1.0;
    double transport_coefficient = 0.0;
    double reaction_coefficient = 0.0;
};

// Struttura per le condizioni al contorno
struct BCConfig {
    enum Type { DIRICHLET, NEUMANN };
    Type type;
    unsigned int tag;
    std::string function;
};

// Struttura per la configurazione del solver
struct SolverConfig {
    double tolerance;
    int max_iterations;
    std::string method;
};

struct QuadratureCfg {
    std::string type = "order2";
};

// Struttura principale di configurazione
struct Config {
    ProblemConfig problem;
    EquationConfig equation;
    std::vector<BCConfig> boundary_conditions;
    SolverConfig solver;
    int quadrature_order = 2;
    template<int dim>
    static std::unique_ptr<QuadratureRule<dim>> make_quadrature(int order) {
        if (order == 2) return std::make_unique<OrderTwoQuadrature<dim>>();
        if (order == 4) {
            if constexpr (dim == 3) {
                throw std::runtime_error("OrderFourQuadrature<3> is not implemented.");
            }
            return std::make_unique<OrderFourQuadrature<dim>>();
        }
        throw std::runtime_error("Quadrature order not supported");
    }
    QuadratureCfg quadrature;

    // Caricamento / validazione / stampa
    static Config loadFromFile(const std::string& filename);
    bool validate() const;
    void print() const;

    // Factory vari
    Grid1D createGrid1D() const;
    Grid2D createGrid2D() const;
    Grid3D createGrid3D() const;

    template<unsigned int dim> Function<dim,1> createForcingFunction() const;
    template<unsigned int dim> Function<dim,1> createDiffusionFunction() const;
    template<unsigned int dim> Function<dim,1> createReactionFunction() const;
    Function<1,1> createTransportFunction1D() const;
    Function<2,2> createTransportFunction2D() const;
    Function<3,3> createTransportFunction3D() const;

    template<unsigned int dim> BoundaryConditions<dim,1> createBoundaryConditions() const;

    template<unsigned int dim>
    std::unique_ptr<QuadratureRule<dim>> createQuadrature() const;
};

#endif // CONFIG_HPP
