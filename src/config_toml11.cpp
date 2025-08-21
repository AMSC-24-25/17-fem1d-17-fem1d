#include "config.hpp"
#include <toml.hpp>  // toml11 library
#include <iostream>
#include <fstream>

Config Config::loadFromFile(const std::string& filename) {
    Config config;
    
    try {
        // Parse TOML file with toml11 - una sola riga!
        const auto data = toml::parse(filename);
        
        // Parse problem section
        const auto& problem = toml::find(data, "problem");
        config.problem.dimension = toml::find_or(problem, "dimension", 2);
        config.problem.mesh_file = toml::find_or(problem, "mesh_file", std::string("mesh/default.msh"));
        config.problem.output_file = toml::find_or(problem, "output_file", std::string("output/solution"));
        config.problem.grid_size = toml::find_or(problem, "grid_size", 100);
        
        // Parse equation section
        const auto& equation = toml::find(data, "equation");
        config.equation.diffusion_coefficient = toml::find_or(equation, "diffusion_coefficient", 1.0);
        config.equation.transport_coefficient = toml::find_or(equation, "transport_coefficient", 0.0);
        config.equation.reaction_coefficient = toml::find_or(equation, "reaction_coefficient", 0.0);
        config.equation.forcing_function = toml::find_or(equation, "forcing_function", std::string("0.0"));
        
        // Parse solver section
        const auto& solver = toml::find(data, "solver");
        config.solver.tolerance = toml::find_or(solver, "tolerance", 1e-12);
        config.solver.max_iterations = toml::find_or(solver, "max_iterations", 1000);
        config.solver.method = toml::find_or(solver, "method", std::string("direct"));
        
        // Parse boundary conditions (arrays) - molto più semplice!
        if (data.contains("boundary_condition")) {
            const auto& bcs = toml::find(data, "boundary_condition").as_array();
            for (const auto& bc_toml : bcs) {
                BCConfig bc;
                bc.tag = toml::find<int>(bc_toml, "tag");
                bc.function = toml::find<std::string>(bc_toml, "function");
                
                std::string type_str = toml::find_or(bc_toml, "type", std::string("dirichlet"));
                bc.type = (type_str == "neumann") ? BCConfig::NEUMANN : BCConfig::DIRICHLET;
                
                config.boundary_conditions.push_back(bc);
            }
        }
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing TOML file '" + filename + "': " + e.what());
    }
    
    if (!config.validate()) {
        throw std::runtime_error("Invalid configuration in file: " + filename);
    }
    
    return config;
}
