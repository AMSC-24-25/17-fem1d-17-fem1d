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
    std::string transport_function_x;
    std::string transport_function_y;
    std::string transport_function_z;
    std::string reaction_function;
    std::string forcing_function;
    
    // Scalar coefficients for backward compatibility
    double diffusion_coefficient = 1.0;
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

struct QuadratureCfg {
    std::string type = "order2";
};

// Main configuration structure
struct Config {
    ProblemConfig problem;
    EquationConfig equation;
    std::vector<BCConfig> boundary_conditions;
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

using symbol_table_t = exprtk::symbol_table<double>;
using expression_t = exprtk::expression<double>;
using parser_t = exprtk::parser<double>;

// Unified thread-safe expression pool for both simple and time-dependent functions
struct ThreadExpressionPool {
    std::vector<double> x_vals, y_vals, z_vals;
    std::vector<double> time_vals;  // Only used for time-dependent expressions
    std::vector<symbol_table_t> symbol_tables;
    std::vector<expression_t> expressions;
    bool is_time_dependent;
    
    ThreadExpressionPool(const std::string& expr_string, unsigned int dimension, int num_threads, bool time_dependent = false);
};

#endif // CONFIG_HPP
