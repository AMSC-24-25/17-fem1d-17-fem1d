/**
 * @file config.hpp
 * @brief Configuration system for FEM problems using TOML files
 * 
 * Provides a configuration framework that parses TOML files
 * to define PDE problems, boundary conditions, and solver parameters.
 * Supports both steady-state and time-dependent problems with mathematical
 * expression parsing using exprtk for coefficient functions.
 */
#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "grid1D.hpp"
#include "grid.hpp"
#include "function.hpp"
#include "boundary_conditions.hpp"
#include "boundary_conditions_td.hpp"
#include "toml.hpp"  // toml11 library (local)
#include "exprtk.hpp"  // exprtk library (local)
#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

// Problem configuration: defines mesh, dimension, and problem type
struct ProblemConfig {
    unsigned int dimension;      // Spatial dimension (1D, 2D, 3D)
    std::string mesh_file;       // Path to mesh file (.msh format)
    double grid1d_start;         // 1D grid start coordinate
    double grid1d_end;           // 1D grid end coordinate
    unsigned int grid1d_size;    // Number of 1D grid points
    std::string output_file;     // Output file prefix
    bool time_dependent = false; // Whether problem is time-dependent
};

// Equation configuration: PDE coefficients as mathematical expressions
struct EquationConfig {
    // Function expressions (e.g., "sin(x)*cos(y)", "1.0 + 0.1*x")
    std::string diffusion_function;     // Diffusion coefficient μ(x,y,z)
    std::string transport_function_x;   // Transport coefficient b_x(x,y,z)
    std::string transport_function_y;   // Transport coefficient b_y(x,y,z)
    std::string transport_function_z;   // Transport coefficient b_z(x,y,z)
    std::string reaction_function;      // Reaction coefficient r(x,y,z)
    std::string forcing_function;       // Forcing term f(x,y,z)
};

// Boundary condition configuration
struct BCConfig {
    enum Type { DIRICHLET, NEUMANN };
    Type type;                         // Condition type
    unsigned int tag;                  // Boundary tag/marker
    std::string function;              // BC function expression
    std::string time_function = "";    // Optional time-dependent function f(x,y,z,t)
};

// Time-dependent problem configuration
struct TimeDependentConfig {
    double final_time = 1.0;               // Simulation end time
    double time_step = 0.01;               // Time step size
    double theta = 0.5;                    // Theta-method parameter: 0=Explicit, 0.5=Crank-Nicolson, 1=Implicit
    std::string initial_condition = "0.0"; // Initial condition u₀(x,y,z)
    std::string forcing_function_td = "";  // Time-dependent forcing f(x,y,z,t)
};

// Quadrature rule configuration
struct QuadratureCfg {
    std::string type = "order2";  // "order2" or "order4"
};

/**
 * @brief Main configuration container for FEM problems
 * 
 * Aggregates all configuration settings and provides factory methods
 * to create configured objects (grids, functions, boundary conditions).
 * Supports TOML file parsing with mathematical expression evaluation.
 */
struct Config {
    ProblemConfig problem;
    EquationConfig equation;
    std::vector<BCConfig> boundary_conditions;
    TimeDependentConfig time_dependent;
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

    // Configuration loading and validation
    static Config loadFromFile(const std::string& filename);
    bool validate() const;
    void print() const;

    // Factory methods for creating configured objects
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

// Expression parsing type aliases
using symbol_table_t = exprtk::symbol_table<double>;
using expression_t = exprtk::expression<double>;
using parser_t = exprtk::parser<double>;

/**
 * @brief Thread-safe expression pool for parallel mathematical expression evaluation
 * 
 * Maintains separate expression instances per thread to avoid race conditions
 * when evaluating mathematical expressions in parallel OpenMP regions.
 */
struct ThreadExpressionPool {
    std::vector<double> x_vals, y_vals, z_vals;        // Spatial coordinates per thread
    std::vector<double> time_vals;                     // Time values (time-dependent only)
    std::vector<symbol_table_t> symbol_tables;         // Symbol tables per thread
    std::vector<expression_t> expressions;             // Compiled expressions per thread
    bool is_time_dependent;                            // Whether expressions include time
    
    ThreadExpressionPool(const std::string& expr_string, unsigned int dimension, int num_threads, bool time_dependent = false);
};

#endif // CONFIG_HPP
