#include "config.hpp"
#include "fem.hpp"
#include "fem_td.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>

// Forward declarations
void solveSteadyStateProblem(const Config& config);
void solveTimeDependentProblem(const Config& config);

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.toml>" << std::endl;
        return 1;
    }
    
    try {
        Config config = Config::loadFromFile(argv[1]);
        std::filesystem::create_directories("../output");
        
        // Check if problem is time-dependent
        if (config.problem.time_dependent) {
            std::cout << "=== TIME-DEPENDENT PROBLEM ===" << std::endl;
            solveTimeDependentProblem(config);
        } else {
            std::cout << "=== STEADY-STATE PROBLEM ===" << std::endl;
            solveSteadyStateProblem(config);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

void solveSteadyStateProblem(const Config& config) {
    if (config.problem.dimension == 1) {
        Grid1D grid = config.createGrid1D();
        Function<1,1> forcing = config.createForcingFunction<1>();
        Function<1,1> diffusion = config.createDiffusionFunction<1>();
        Function<1,1> transport = config.createTransportFunction1D();
        Function<1,1> reaction = config.createReactionFunction<1>();
        BoundaryConditions<1,1> bc = config.createBoundaryConditions<1>();
        std::unique_ptr<QuadratureRule<1>> quadrature = config.createQuadrature<1>();
        
        Fem<1> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
        fem.assemble();
        fem.solve();
        fem.outputCsv("../" + config.problem.output_file + ".csv");
        fem.outputVtu("../" + config.problem.output_file + ".vtu");
        
    } else if (config.problem.dimension == 2) {
        Grid2D grid = config.createGrid2D();
        Function<2,1> forcing = config.createForcingFunction<2>();
        Function<2,1> diffusion = config.createDiffusionFunction<2>();
        Function<2,2> transport = config.createTransportFunction2D();
        Function<2,1> reaction = config.createReactionFunction<2>();
        BoundaryConditions<2,1> bc = config.createBoundaryConditions<2>();
        std::unique_ptr<QuadratureRule<2>> quadrature = config.createQuadrature<2>();
        
        Fem<2> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
        fem.assemble();
        fem.solve();
        fem.outputVtu("../" + config.problem.output_file + ".vtu");
        fem.outputCsv("../" + config.problem.output_file + ".csv");
        
    } else if (config.problem.dimension == 3) {
        Grid3D grid = config.createGrid3D();
        Function<3,1> forcing = config.createForcingFunction<3>();
        Function<3,1> diffusion = config.createDiffusionFunction<3>();
        Function<3,3> transport = config.createTransportFunction3D();
        Function<3,1> reaction = config.createReactionFunction<3>();
        BoundaryConditions<3,1> bc = config.createBoundaryConditions<3>();
        std::unique_ptr<QuadratureRule<3>> quadrature = config.createQuadrature<3>();
        
        Fem<3> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
        fem.assemble();
        fem.solve();
        fem.outputVtu("../" + config.problem.output_file + ".vtu");
        fem.outputCsv("../" + config.problem.output_file + ".csv");
        
    } else {
        std::cerr << "Error: Unsupported dimension " << config.problem.dimension << std::endl;
    }
}

void solveTimeDependentProblem(const Config& config) {
    if (config.problem.dimension == 1) {
        Grid1D grid = config.createGrid1D();
        Function<1,1> diffusion = config.createDiffusionFunction<1>();
        Function<1,1> transport = config.createTransportFunction1D();
        Function<1,1> reaction = config.createReactionFunction<1>();
        BoundaryConditions_td<1,1> bc = config.createBoundaryConditionsTD<1>();
        std::unique_ptr<QuadratureRule<1>> quadrature = config.createQuadrature<1>();
        
        FemTD<1> femtd(grid, diffusion, transport, reaction, bc, *quadrature);
        
        // Set time-dependent forcing (simplified)
        femtd.set_forcing([](const Point<1>& p, double t) -> double {
            return std::sin(2*M_PI*p[0]) * std::cos(t);
        });
        
        // Set initial condition
        Function<1,1> initial([](const Point<1>& p) -> double { return 0.0; });
        femtd.set_initial_condition(initial);
        
        // Run simulation
        femtd.run(config.time_dependent.final_time, 
                  config.time_dependent.time_step, 
                  config.time_dependent.theta,
                  "../" + config.problem.output_file);
                  
    } else if (config.problem.dimension == 2) {
        Grid2D grid = config.createGrid2D();
        Function<2,1> diffusion = config.createDiffusionFunction<2>();
        Function<2,2> transport = config.createTransportFunction2D();
        Function<2,1> reaction = config.createReactionFunction<2>();
        BoundaryConditions_td<2,1> bc = config.createBoundaryConditionsTD<2>();
        std::unique_ptr<QuadratureRule<2>> quadrature = config.createQuadrature<2>();
        
        FemTD<2> femtd(grid, diffusion, transport, reaction, bc, *quadrature);
        
        // Set time-dependent forcing (simplified)
        femtd.set_forcing([](const Point<2>& p, double t) -> double {
            return std::sin(2*M_PI*p[0]*p[1]) * std::exp(-t);
        });
        
        // Set initial condition
        Function<2,1> initial([](const Point<2>& p) -> double {
            return std::exp(-((p[0]-0.5)*(p[0]-0.5) + (p[1]-0.5)*(p[1]-0.5))/0.1);
        });
        femtd.set_initial_condition(initial);
        
        // Pre-assemble time-invariant matrices
        femtd.assemble_time_invariant();
        
        // Run simulation
        femtd.run(config.time_dependent.final_time, 
                  config.time_dependent.time_step, 
                  config.time_dependent.theta,
                  "../" + config.problem.output_file);
                  
    } else if (config.problem.dimension == 3) {
        Grid3D grid = config.createGrid3D();
        Function<3,1> diffusion = config.createDiffusionFunction<3>();
        Function<3,3> transport = config.createTransportFunction3D();
        Function<3,1> reaction = config.createReactionFunction<3>();
        BoundaryConditions_td<3,1> bc = config.createBoundaryConditionsTD<3>();
        std::unique_ptr<QuadratureRule<3>> quadrature = config.createQuadrature<3>();
        
        FemTD<3> femtd(grid, diffusion, transport, reaction, bc, *quadrature);
        
        // Set time-dependent forcing (simplified)
        femtd.set_forcing([](const Point<3>& p, double t) -> double {
            return std::sin(2*M_PI*p[0]*p[1]*p[2]) + 4*M_PI*M_PI*t*std::sin(2*M_PI*p[0]*p[1]*p[2]);
        });
        
        // Set initial condition
        Function<3,1> initial([](const Point<3>& p) -> double { return 0.0; });
        femtd.set_initial_condition(initial);
        
        // Pre-assemble time-invariant matrices
        femtd.assemble_time_invariant();
        
        // Run simulation
        femtd.run(config.time_dependent.final_time, 
                  config.time_dependent.time_step, 
                  config.time_dependent.theta,
                  "../" + config.problem.output_file);
                  
    } else {
        std::cerr << "Error: Unsupported dimension " << config.problem.dimension << std::endl;
    }
}
