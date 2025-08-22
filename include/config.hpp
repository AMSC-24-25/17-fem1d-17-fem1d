#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>
#include "grid1D.hpp"
#include "grid.hpp"
#include "function.hpp"
#include "boundary_conditions.hpp"
#include "boundary_conditions_td.hpp"
#include <memory>

// Structure for problem configuration
struct ProblemConfig {
    unsigned int dimension;
    std::string mesh_file;
    std::string output_file;
    int grid_size;  // For 1D uniform grids
    bool time_dependent = false;  // NEW: flag for time-dependent problems
};

// Structure for equation configuration
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

// Structure for boundary condition configuration
struct BCConfig {
    enum Type { DIRICHLET, NEUMANN };
    Type type;
    unsigned int tag;
    std::string function;
    std::string time_function = "";  // NEW: optional time-dependent function "f(x,y,z,t)"
};

// NEW: Structure for time-dependent configuration
struct TimeDependentConfig {
    double final_time = 1.0;
    double time_step = 0.01;
    double theta = 0.5;  // 0=Explicit, 0.5=Crank-Nicolson, 1=Implicit
    std::string initial_condition = "0.0";
    std::string forcing_function_td = "";  // f(x,y,z,t) - time dependent forcing
};

// Structure for solver configuration
struct SolverConfig {
    double tolerance;
    int max_iterations;
    std::string method;
};

struct QuadratureCfg {
    std::string type = "order2";
};

// Main configuration structure
struct Config {
    ProblemConfig problem;
    EquationConfig equation;
    std::vector<BCConfig> boundary_conditions;
    SolverConfig solver;
    TimeDependentConfig time_dependent;  // NEW: time-dependent settings
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

    // Loading / validation / printing
    static Config loadFromFile(const std::string& filename);
    bool validate() const;
    void print() const;

    // Various factories
    template<unsigned int dim>
    Grid<dim> createGrid() const;

    template<unsigned int dim> Function<dim,1> createForcingFunction() const;
    template<unsigned int dim>
    std::function<double(const Point<dim>&, double)> createForcingFunction_td() const;

    template<unsigned int dim> Function<dim,1> createDiffusionFunction() const;
    template<unsigned int dim> Function<dim,1> createReactionFunction() const;
    template<unsigned int dim> Function<dim,dim> createTransportFunction() const;

    template<unsigned int dim> BoundaryConditions<dim,1> createBoundaryConditions() const;
    template<unsigned int dim> BoundaryConditions_td<dim,1> createBoundaryConditionsTD() const; 

    template<unsigned int dim> std::unique_ptr<QuadratureRule<dim>> createQuadrature() const;

    template<unsigned int dim> Function<dim,1> createInitialConditionFunction() const;
};

#endif // CONFIG_HPP
