#include "config.hpp"
#include "toml.hpp"  //         config.equation.forcing_function = toml::find_or(equation, "forcing_function", std::string("0.0"));
        
        // Parse boundary conditions (arrays) - much more elegant!rary (local)
#include "exprtk.hpp"  // exprtk library (local)
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

Config Config::loadFromFile(const std::string& filename) {
    Config config;
    
    try {
    // Parse TOML file with toml11 v4 - much simpler!
        const toml::value data = toml::parse(filename);
        
        // Parse problem section
        const toml::value& problem = toml::find(data, "problem");
        config.problem.dimension = toml::find_or(problem, "dimension", 2);
        config.problem.mesh_file = toml::find_or(problem, "mesh_file", std::string("mesh/default.msh"));
        config.problem.output_file = toml::find_or(problem, "output_file", std::string("output/solution"));
        config.problem.grid_size = toml::find_or(problem, "grid_size", 100);
        config.problem.time_dependent = toml::find_or(problem, "time_dependent", false);  // NEW
        
    // Parse equation section - unified approach: everything is a function
        const toml::value& equation = toml::find(data, "equation");
        
    // Read as functions first, fallback to coefficient values
        config.equation.diffusion_function = toml::find_or(equation, "diffusion_function", 
            toml::find_or(equation, "diffusion_coefficient", std::string("1.0")));
        config.equation.transport_function_x = toml::find_or(equation, "transport_function_x", 
            toml::find_or(equation, "transport_coefficient_x", std::string("0.0")));
        config.equation.transport_function_y = toml::find_or(equation, "transport_function_y", 
            toml::find_or(equation, "transport_coefficient_y", std::string("0.0")));
        config.equation.transport_function_z = toml::find_or(equation, "transport_function_z", 
            toml::find_or(equation, "transport_coefficient_z", std::string("0.0")));
        config.equation.reaction_function = toml::find_or(equation, "reaction_function", 
            toml::find_or(equation, "reaction_coefficient", std::string("0.0")));
        config.equation.forcing_function = toml::find_or(equation, "forcing_function", std::string("0.0"));
        
        // Parse boundary conditions (arrays) - much more elegant!
        if (data.contains("boundary_conditions")) {
            const toml::array& bcs = toml::find(data, "boundary_conditions").as_array();
            for (const toml::value& bc_toml : bcs) {
                BCConfig bc;
                bc.tag = toml::find<int>(bc_toml, "tag");
                bc.function = toml::find_or<std::string>(bc_toml, "function", std::string("0.0"));
                bc.time_function = toml::find_or(bc_toml, "time_function", std::string("")); 
                
                std::string type_str = toml::find_or(bc_toml, "type", std::string("dirichlet"));
                bc.type = (type_str == "neumann") ? BCConfig::NEUMANN : BCConfig::DIRICHLET;
                
                config.boundary_conditions.push_back(bc);
            }
        }

        // NEW: Parse time-dependent section
        if (data.contains("time_dependent")) {
            const toml::value& td = toml::find(data, "time_dependent");
            config.time_dependent.final_time = toml::find_or(td, "final_time", 1.0);
            config.time_dependent.time_step = toml::find_or(td, "time_step", 0.01);
            config.time_dependent.theta = toml::find_or(td, "theta", 0.5);
            config.time_dependent.initial_condition = toml::find_or(td, "initial_condition", std::string("0.0"));
            config.time_dependent.forcing_function_td = toml::find_or(td, "forcing_function_td", std::string(""));
        }
        if (data.contains("quadrature")) {
            const toml::value& quad = toml::find(data, "quadrature");
            std::string qtype = toml::find_or(quad, "type", std::string("order2"));
            std::transform(qtype.begin(), qtype.end(), qtype.begin(), ::tolower);
            config.quadrature.type = qtype;
        } else {
            config.quadrature.type = "order2"; // default
        }
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing TOML file '" + filename + "': " + e.what());
    }
    
    if (!config.validate()) {
        throw std::runtime_error("Invalid configuration in file: " + filename);
    }
    
    return config;
}

bool Config::validate() const {
    // Validate dimension
    if (problem.dimension < 1 || problem.dimension > 3) {
        std::cerr << "Error: dimension must be 1, 2, or 3" << std::endl;
        return false;
    }
    
    // Validate mesh file
    if (problem.mesh_file.empty()) {
        std::cerr << "Error: mesh_file cannot be empty" << std::endl;
        return false;
    }
    
    // Validate boundary conditions
    if (boundary_conditions.empty()) {
        std::cerr << "Warning: no boundary conditions specified" << std::endl;
    }

    if (!(quadrature.type == "order2" || quadrature.type == "order4")) {
        throw std::runtime_error("quadrature.type must be 'order2' or 'order4'. Got: " + quadrature.type);
    }
    
    return true;
}

void Config::print() const {
    std::cout << "=== Configuration ===" << std::endl;
    std::cout << "Problem:" << std::endl;
    std::cout << "  Dimension: " << problem.dimension << std::endl;
    std::cout << "  Mesh file: " << problem.mesh_file << std::endl;
    std::cout << "  Output file: " << problem.output_file << std::endl;
    std::cout << "  Time-dependent: " << (problem.time_dependent ? "Yes" : "No") << std::endl;  // NEW
    
    std::cout << "Equation:" << std::endl;
    std::cout << "  Diffusion function: " << equation.diffusion_function << std::endl;
    std::cout << "  Transport function (x): " << equation.transport_function_x << std::endl;
    std::cout << "  Transport function (y): " << equation.transport_function_y << std::endl;
    std::cout << "  Transport function (z): " << equation.transport_function_z << std::endl;
    std::cout << "  Reaction function: " << equation.reaction_function << std::endl;
    std::cout << "  Forcing function: " << equation.forcing_function << std::endl;
    
    std::cout << "Boundary Conditions:" << std::endl;
    for (size_t i = 0; i < boundary_conditions.size(); ++i) {
        const BCConfig& bc = boundary_conditions[i];
        std::cout << "  BC " << i+1 << ": ";
        std::cout << (bc.type == BCConfig::DIRICHLET ? "Dirichlet" : "Neumann");
        std::cout << " on tag " << bc.tag;
        std::cout << ", function: " << bc.function;
        if (!bc.time_function.empty()) {
            std::cout << ", time_function: " << bc.time_function;  // NEW
        }
        std::cout << std::endl;
    }
    
    // NEW: Print time-dependent config if applicable
    if (problem.time_dependent) {
        std::cout << "Time-Dependent Settings:" << std::endl;
        std::cout << "  Final time: " << time_dependent.final_time << std::endl;
        std::cout << "  Time step: " << time_dependent.time_step << std::endl;
        std::cout << "  Theta: " << time_dependent.theta << std::endl;
        std::cout << "  Initial condition: " << time_dependent.initial_condition << std::endl;
        if (!time_dependent.forcing_function_td.empty()) {
            std::cout << "  TD Forcing function: " << time_dependent.forcing_function_td << std::endl;
        }
    }
    
    std::cout << "===================" << std::endl;
}

ThreadExpressionPool::ThreadExpressionPool(const std::string& expr_string, unsigned int dimension, int num_threads, bool time_dependent = false) 
  : x_vals(num_threads, 0.0),
    y_vals(num_threads, 0.0),
    z_vals(num_threads, 0.0),
    is_time_dependent(time_dependent) {
        
    if (is_time_dependent) {
        time_vals.resize(num_threads, 0.0);
    }
    symbol_tables.resize(num_threads);
    expressions.resize(num_threads);
    
    // Initialize each thread's expression
    for (int i = 0; i < num_threads; ++i) {
        symbol_tables[i].add_variable("x", x_vals[i]);
        if (dimension >= 2) symbol_tables[i].add_variable("y", y_vals[i]);
        if (dimension >= 3) symbol_tables[i].add_variable("z", z_vals[i]);
        
        // Add time variables only for time-dependent expressions
        if (is_time_dependent) {
            symbol_tables[i].add_variable("t", time_vals[i]);
            symbol_tables[i].add_variable("time", time_vals[i]);
        }
        
        symbol_tables[i].add_constant("pi", M_PI);
        symbol_tables[i].add_constant("e", M_E);
        
        expressions[i].register_symbol_table(symbol_tables[i]);
        
        parser_t parser;
        if (!parser.compile(expr_string, expressions[i])) {
            throw std::runtime_error("Failed to compile expression: " + expr_string);
        }
    }
}

// Enhanced function parsing with exprtk - supports any mathematical expression!
template<unsigned int dim>
Function<dim,1> parseSimpleFunction(const std::string& expression) {
    // Optimization for simple constant functions
    if (expression == "0" || expression == "0.0") {
        return Function<dim,1>([](Point<dim> p) -> double { return 0.0; });
    }
    if (expression == "1" || expression == "1.0") {
        return Function<dim,1>([](Point<dim> p) -> double { return 1.0; });
    }

    // Thread pool approach: create one expression per potential thread
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
#else
    const int max_threads = 1;
#endif

    std::shared_ptr<ThreadExpressionPool> pool = std::make_shared<ThreadExpressionPool>(
                expression, dim, max_threads, false // false = not time-dependent
    );

    // Return thread-safe function using thread pool
    return Function<dim,1>([pool](Point<dim> p) -> double {
        #ifdef _OPENMP
            const int thread_id = omp_get_thread_num();
        #else
            const int thread_id = 0;
        #endif
        // Update variables for this thread
        pool->x_vals[thread_id] = (dim >= 1) ? p[0] : 0.0;
        pool->y_vals[thread_id] = (dim >= 2) ? p[1] : 0.0;
        pool->z_vals[thread_id] = (dim >= 3) ? p[2] : 0.0;
        
        return pool->expressions[thread_id].value();
    });
}

// NEW: Parse time-dependent function f(x,y,z,t)
template<unsigned int dim>
std::function<double(const Point<dim>&, double)> parseTimeDependentFunction(const std::string& expression) {
    if (expression.empty() || expression == "0" || expression == "0.0") {
        return [](const Point<dim>&, double) -> double { return 0.0; };
    }
    if (expression == "1" || expression == "1.0") {
        return [](const Point<dim>&, double) -> double { return 1.0; };
    }

    // Thread pool approach: create one expression per potential thread
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
#else
    const int max_threads = 1;
#endif

    std::shared_ptr<ThreadExpressionPool> pool = std::make_shared<ThreadExpressionPool>(
                expression, dim, max_threads, true // true = time-dependent
    ); 

    // Return thread-safe function using thread pool
    return [pool](const Point<dim>& p, double t) -> double {
        #ifdef _OPENMP
            const int thread_id = omp_get_thread_num();
        #else
            const int thread_id = 0;
        #endif
        
        // Update variables for this thread
        pool->x_vals[thread_id] = (dim >= 1) ? p[0] : 0.0;
        pool->y_vals[thread_id] = (dim >= 2) ? p[1] : 0.0;
        pool->z_vals[thread_id] = (dim >= 3) ? p[2] : 0.0;
        pool->time_vals[thread_id] = t;  // Set time value
        
        return pool->expressions[thread_id].value();
    };
}

// Config factory methods implementation
template<>
Grid<1> Config::createGrid<1>() const {
    if (!problem.mesh_file.empty() && problem.mesh_file != "mesh/default.msh") {
        // If a specific mesh file is provided, try to load it
        std::cerr << "Warning: 1D mesh file loading not implemented, using uniform grid" << std::endl;
    }
    // Create uniform grid from 0 to 1 with grid_size points
    return Grid1D(0.0, 1.0, problem.grid_size);
}

template<>
Grid<2> Config::createGrid<2>() const {
    Grid<2> grid;
    grid.parseFromMsh(problem.mesh_file);
    return grid;
}

template<>
Grid<3> Config::createGrid<3>() const {
    Grid<3> grid;
    grid.parseFromMsh(problem.mesh_file);
    return grid;
}

template<unsigned int dim>
Function<dim,1> Config::createForcingFunction() const {
    return parseSimpleFunction<dim>(equation.forcing_function);
}
template<unsigned int dim>
std::function<double(const Point<dim>&, double)> Config::createForcingFunction_td() const {
    return parseTimeDependentFunction<dim>(time_dependent.forcing_function_td);
}
template<unsigned int dim>
Function<dim,1> Config::createInitialConditionFunction() const {
    return parseSimpleFunction<dim>(time_dependent.initial_condition);
}

template<unsigned int dim>
Function<dim,1> Config::createDiffusionFunction() const {
    return parseSimpleFunction<dim>(equation.diffusion_function);
}

template<unsigned int dim>
Function<dim,1> Config::createReactionFunction() const {
    return parseSimpleFunction<dim>(equation.reaction_function);
}

template<>
Function<1,1> Config::createTransportFunction<1>() const {
    return parseSimpleFunction<1>(equation.transport_function_x);
}

template<>
Function<2,2> Config::createTransportFunction<2>() const {
    Function<2,1> xFunc = parseSimpleFunction<2>(equation.transport_function_x);
    Function<2,1> yFunc = parseSimpleFunction<2>(equation.transport_function_y);
    return Function<2,2>([xFunc, yFunc](Point<2> p) { 
        double xCoeff = xFunc(p);
        double yCoeff = yFunc(p);
        return Point<2>(xCoeff, yCoeff);
    });
}

template<>
Function<3,3> Config::createTransportFunction<3>() const {
    Function<3,1> xFunc = parseSimpleFunction<3>(equation.transport_function_x);
    Function<3,1> yFunc = parseSimpleFunction<3>(equation.transport_function_y);
    Function<3,1> zFunc = parseSimpleFunction<3>(equation.transport_function_z);
    return Function<3,3>([xFunc, yFunc, zFunc](Point<3> p) { 
        double xCoeff = xFunc(p);
        double yCoeff = yFunc(p);
        double zCoeff = zFunc(p);
        return Point<3>(xCoeff, yCoeff, zCoeff);
    });
}

template<unsigned int dim>
BoundaryConditions<dim,1> Config::createBoundaryConditions() const {
    BoundaryConditions<dim,1> bc;

    for (const BCConfig& bcConfig : boundary_conditions) {
        Function<dim,1> func = parseSimpleFunction<dim>(bcConfig.function);
        
        if (bcConfig.type == BCConfig::DIRICHLET) {
            bc.addDirichlet(bcConfig.tag, func);
        } else if (bcConfig.type == BCConfig::NEUMANN) {
            bc.addNeumann(bcConfig.tag, func);
        }
    }
    
    return bc;
}

// NEW: Create time-dependent boundary conditions
template<unsigned int dim>
BoundaryConditions_td<dim,1> Config::createBoundaryConditionsTD() const {
    BoundaryConditions_td<dim,1> bc;

    for (const BCConfig& bcConfig : boundary_conditions) {
        // Use time_function if available, otherwise fall back to regular function (time-independent)
        std::string function_expr = bcConfig.time_function.empty() ? bcConfig.function : bcConfig.time_function;
        
        if (bcConfig.time_function.empty()) {
            // Time-independent: convert regular function to time-dependent
            Function<dim,1> static_func = parseSimpleFunction<dim>(bcConfig.function);
            fun_td<dim,1> td_func = [static_func](const Point<dim>& p, double t) -> double {
                return static_func.value(p);
            };
            
            if (bcConfig.type == BCConfig::DIRICHLET) {
                bc.addDirichlet(bcConfig.tag, td_func);
            } else if (bcConfig.type == BCConfig::NEUMANN) {
                bc.addNeumann(bcConfig.tag, td_func);
            }
        } else {
            // Time-dependent function
            fun_td<dim,1> td_func = parseTimeDependentFunction<dim>(bcConfig.time_function);
            
            if (bcConfig.type == BCConfig::DIRICHLET) {
                bc.addDirichlet(bcConfig.tag, td_func);
            } else if (bcConfig.type == BCConfig::NEUMANN) {
                bc.addNeumann(bcConfig.tag, td_func);
            }
        }
    }
    
    return bc;
}

template<>
std::unique_ptr<QuadratureRule<1>> Config::createQuadrature<1>() const {
    if (quadrature.type == "order4")
        return std::make_unique<OrderFourQuadrature<1>>();
    // safe default: order2
    return std::make_unique<OrderTwoQuadrature<1>>();
}

template<>
std::unique_ptr<QuadratureRule<2>> Config::createQuadrature<2>() const {
    if (quadrature.type == "order4")
        return std::make_unique<OrderFourQuadrature<2>>();
    return std::make_unique<OrderTwoQuadrature<2>>();
}

template<>
std::unique_ptr<QuadratureRule<3>> Config::createQuadrature<3>() const {
    if (quadrature.type == "order4") {
        throw std::runtime_error("OrderFourQuadrature<3> is not implemented.");
    }
    return std::make_unique<OrderTwoQuadrature<3>>();
}


#define INSTANTIATE_CONFIG_FUNCS(DIM) \
    template Grid<DIM> Config::createGrid<DIM>() const; \
    template Function<DIM,1> Config::createForcingFunction<DIM>() const; \
    template std::function<double(const Point<DIM>&, double)> Config::createForcingFunction_td<DIM>() const; \
    template Function<DIM,1> Config::createDiffusionFunction<DIM>() const; \
    template Function<DIM,1> Config::createReactionFunction<DIM>() const; \
    template Function<DIM,1> Config::createInitialConditionFunction<DIM>() const; \
    template Function<DIM,DIM> Config::createTransportFunction<DIM>() const; \
    template BoundaryConditions<DIM,1> Config::createBoundaryConditions<DIM>() const; \
    template BoundaryConditions_td<DIM,1> Config::createBoundaryConditionsTD<DIM>() const; \
    template std::unique_ptr<QuadratureRule<DIM>> Config::createQuadrature<DIM>() const;

INSTANTIATE_CONFIG_FUNCS(1)
INSTANTIATE_CONFIG_FUNCS(2)
INSTANTIATE_CONFIG_FUNCS(3)