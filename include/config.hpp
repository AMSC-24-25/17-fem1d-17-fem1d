#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>
#include "grid1D.hpp"
#include "grid2D.hpp"
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

// Struttura principale di configurazione
struct Config {
    ProblemConfig problem;
    EquationConfig equation;
    std::vector<BCConfig> boundary_conditions;
    SolverConfig solver;
    
    // Metodo statico per caricare da file TOML
    static Config loadFromFile(const std::string& filename);
    
    // Metodo per validare la configurazione
    bool validate() const;
    
    // Metodo per stampare la configurazione (debug)
    void print() const;
    
    // Factory methods for creating FEM objects
    Grid1D createGrid1D() const;
    Grid2D createGrid2D() const;
    
    template<unsigned int dim>
    Function<dim,1> createForcingFunction() const;
    
    template<unsigned int dim>
    Function<dim,1> createDiffusionFunction() const;
    
    template<unsigned int dim>
    Function<dim,1> createReactionFunction() const;
    
    Function<1,1> createTransportFunction1D() const;
    Function<2,2> createTransportFunction2D() const;
    
    template<unsigned int dim>
    BoundaryConditions<dim,1> createBoundaryConditions() const;
};

#endif // CONFIG_HPP
