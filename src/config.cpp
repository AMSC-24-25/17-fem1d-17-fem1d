#include "config.hpp"
#include "toml.hpp"  // toml11 library (local)
#include "exprtk.hpp"  // exprtk library (local)
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Config Config::loadFromFile(const std::string& filename) {
    Config config;
    
    try {
    // Parse TOML file with toml11 v4 - much simpler!
        //from auto to toml::value
        const toml::value data = toml::parse(filename);
        
        // Parse problem section
        //from const auto& to const toml::value&
        const toml::value& problem = toml::find(data, "problem");
        config.problem.dimension = toml::find_or(problem, "dimension", 2);
        config.problem.mesh_file = toml::find_or(problem, "mesh_file", std::string("mesh/default.msh"));
        config.problem.output_file = toml::find_or(problem, "output_file", std::string("output/solution"));
        config.problem.grid_size = toml::find_or(problem, "grid_size", 100);
        config.problem.time_dependent = toml::find_or(problem, "time_dependent", false);  // NEW
        
    // Parse equation section - unified approach: everything is a function
        //from const auto& to const toml::value&
        const toml::value& equation = toml::find(data, "equation");
        
    // Read as functions first, fallback to coefficient values
        config.equation.diffusion_function = toml::find_or(equation, "diffusion_function", 
            toml::find_or(equation, "diffusion_coefficient", std::string("1.0")));
        config.equation.transport_function = toml::find_or(equation, "transport_function", 
            toml::find_or(equation, "transport_coefficient", std::string("0.0")));
        config.equation.reaction_function = toml::find_or(equation, "reaction_function", 
            toml::find_or(equation, "reaction_coefficient", std::string("0.0")));
        config.equation.forcing_function = toml::find_or(equation, "forcing_function", std::string("0.0"));
        
    // Parse solver section
        //from const auto& to const toml::value&
        const toml::value& solver = toml::find(data, "solver");
        config.solver.tolerance = toml::find_or(solver, "tolerance", 1e-12);
        config.solver.max_iterations = toml::find_or(solver, "max_iterations", 1000);
        config.solver.method = toml::find_or(solver, "method", std::string("direct"));
        
        // Parse boundary conditions (arrays) - much more elegant!
        if (data.contains("boundary_conditions")) {
            //from const auto& to const toml::array&
            const toml::array& bcs = toml::find(data, "boundary_conditions").as_array();
            //from const auto& to const toml::value&
            for (const toml::value& bc_toml : bcs) {
                BCConfig bc;
                bc.tag = toml::find<int>(bc_toml, "tag");
                bc.function = toml::find<std::string>(bc_toml, "function");
                bc.time_function = toml::find_or(bc_toml, "time_function", std::string(""));  // NEW
                
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
        }        if (data.contains("quadrature")) {
            //from const auto& to const toml::value&
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
    
    // Validate solver method
    if (solver.method != "direct" && solver.method != "iterative") {
        std::cerr << "Error: solver method must be 'direct' or 'iterative'" << std::endl;
        return false;
    }
    
    // Validate tolerance
    if (solver.tolerance <= 0) {
        std::cerr << "Error: tolerance must be positive" << std::endl;
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
    std::cout << "  Transport function: " << equation.transport_function << std::endl;
    std::cout << "  Reaction function: " << equation.reaction_function << std::endl;
    std::cout << "  Forcing function: " << equation.forcing_function << std::endl;
    
    std::cout << "Solver:" << std::endl;
    std::cout << "  Tolerance: " << solver.tolerance << std::endl;
    std::cout << "  Max iterations: " << solver.max_iterations << std::endl;
    std::cout << "  Method: " << solver.method << std::endl;
    
    std::cout << "Boundary Conditions:" << std::endl;
    for (size_t i = 0; i < boundary_conditions.size(); ++i) {
        //from const auto& to const BCConfig&
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
    
    // For all other expressions, use exprtk directly
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    
    return Function<dim,1>([expression](Point<dim> p) -> double {
        symbol_table_t symbol_table;
        expression_t expr;
        parser_t parser;
        
        // Register variables based on dimension
        double x = (dim >= 1) ? p[0] : 0.0;
        double y = (dim >= 2) ? p[1] : 0.0;
        double z = (dim >= 3) ? p[2] : 0.0;
        
        symbol_table.add_variable("x", x);
        if (dim >= 2) symbol_table.add_variable("y", y);
        if (dim >= 3) symbol_table.add_variable("z", z);
        
        // Add common constants
        symbol_table.add_constant("pi", M_PI);
        symbol_table.add_constant("e", M_E);
        
        expr.register_symbol_table(symbol_table);
        
        if (parser.compile(expression, expr)) {
            return expr.value();
        } else {
            std::cerr << "Error parsing expression '" << expression << "': " << parser.error() << std::endl;
            return 0.0;
        }
    });
}

// NEW: Parse time-dependent function f(x,y,z,t)
template<unsigned int dim>
std::function<double(const Point<dim>&, double)> parseTimeDependentFunction(const std::string& expression) {
    if (expression.empty() || expression == "0" || expression == "0.0") {
        return [](const Point<dim>& p, double t) -> double { return 0.0; };
    }
    
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    
    return [expression](const Point<dim>& p, double t) -> double {
        symbol_table_t symbol_table;
        expression_t expr;
        parser_t parser;
        
        // Register variables based on dimension
        double x = (dim >= 1) ? p[0] : 0.0;
        double y = (dim >= 2) ? p[1] : 0.0;
        double z = (dim >= 3) ? p[2] : 0.0;
        double time = t;
        
        symbol_table.add_variable("x", x);
        if (dim >= 2) symbol_table.add_variable("y", y);
        if (dim >= 3) symbol_table.add_variable("z", z);
        symbol_table.add_variable("t", time);
        symbol_table.add_variable("time", time);
        
        // Add common constants
        symbol_table.add_constant("pi", M_PI);
        symbol_table.add_constant("e", M_E);
        
        expr.register_symbol_table(symbol_table);
        
        if (parser.compile(expression, expr)) {
            return expr.value();
        } else {
            std::cerr << "Error parsing time-dependent expression '" << expression << "': " << parser.error() << std::endl;
            return 0.0;
        }
    };
}

// Config factory methods implementation
Grid1D Config::createGrid1D() const {
    if (!problem.mesh_file.empty() && problem.mesh_file != "mesh/default.msh") {
        // If a specific mesh file is provided, try to load it
        std::cerr << "Warning: 1D mesh file loading not implemented, using uniform grid" << std::endl;
    }
    // Create uniform grid from 0 to 1 with grid_size points
    return Grid1D(0.0, 1.0, problem.grid_size);
}

Grid2D Config::createGrid2D() const {
    Grid2D grid;
    grid.parseFromMsh(problem.mesh_file);
    return grid;
}

Grid3D Config::createGrid3D() const {
    Grid3D grid;
    grid.parseFromMsh(problem.mesh_file);
    return grid;
}

template<unsigned int dim>
Function<dim,1> Config::createForcingFunction() const {
    return parseSimpleFunction<dim>(equation.forcing_function);
}
template<unsigned int dim>
std::function<double(const Point<dim>&, double)> Config::createForcingFunction_td() const {
    return parseTimeDependentFunction<dim>(equation.forcing_function);
}
template<unsigned int dim>
Function<dim,1> Config::createInitialConditionFunction() const {
    return parseSimpleFunction<dim>(equation.initial_condition);
}

template<unsigned int dim>
Function<dim,1> Config::createDiffusionFunction() const {
    return parseSimpleFunction<dim>(equation.diffusion_function);
}

template<unsigned int dim>
Function<dim,1> Config::createReactionFunction() const {
    return parseSimpleFunction<dim>(equation.reaction_function);
}

Function<1,1> Config::createTransportFunction1D() const {
    return parseSimpleFunction<1>(equation.transport_function);
}

Function<2,2> Config::createTransportFunction2D() const {
    // Parse as scalar and assume transport in x direction
    // from auto to Function<2,1>
    Function<2,1> scalarFunc = parseSimpleFunction<2>(equation.transport_function);
    return Function<2,2>([scalarFunc](Point<2> p) { 
        double coeff = scalarFunc(p);
        return Point<2>(coeff, 0.0);  // transport in x direction
    });
}

Function<3,3> Config::createTransportFunction3D() const {
    // Parse as scalar and assume transport in x direction
    //from auto to Function<3,1>
    Function<3,1> scalarFunc = parseSimpleFunction<3>(equation.transport_function);
    return Function<3,3>([scalarFunc](Point<3> p) { 
        double coeff = scalarFunc(p);
        return Point<3>(coeff, 0.0, 0.0);  // transport in x direction
    });
}

template<unsigned int dim>
BoundaryConditions<dim,1> Config::createBoundaryConditions() const {
    BoundaryConditions<dim,1> bc;

    //from const auto& to const BCConfig&
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
            auto td_func = [static_func](const Point<dim>& p, double t) -> double {
                return static_func.value(p);
            };
            
            if (bcConfig.type == BCConfig::DIRICHLET) {
                bc.addDirichlet(bcConfig.tag, td_func);
            } else if (bcConfig.type == BCConfig::NEUMANN) {
                bc.addNeumann(bcConfig.tag, td_func);
            }
        } else {
            // Time-dependent function
            auto td_func = parseTimeDependentFunction<dim>(bcConfig.time_function);
            
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

// Explicit template instantiations
template Function<1,1> Config::createForcingFunction<1>() const;
template Function<2,1> Config::createForcingFunction<2>() const;
template Function<3,1> Config::createForcingFunction<3>() const;
template Function<1,1> Config::createDiffusionFunction<1>() const;
template Function<2,1> Config::createDiffusionFunction<2>() const;
template Function<3,1> Config::createDiffusionFunction<3>() const;
template Function<1,1> Config::createReactionFunction<1>() const;
template Function<2,1> Config::createReactionFunction<2>() const;
template Function<3,1> Config::createReactionFunction<3>() const;
template BoundaryConditions<1,1> Config::createBoundaryConditions<1>() const;
template BoundaryConditions<2,1> Config::createBoundaryConditions<2>() const;
template BoundaryConditions<3,1> Config::createBoundaryConditions<3>() const;
template BoundaryConditions_td<1,1> Config::createBoundaryConditionsTD<1>() const;
template BoundaryConditions_td<2,1> Config::createBoundaryConditionsTD<2>() const;
template BoundaryConditions_td<3,1> Config::createBoundaryConditionsTD<3>() const;