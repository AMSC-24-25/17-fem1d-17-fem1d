#include "config.hpp"
#include "fem.hpp"
#include "fem_td.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>

// Forward declarations
template<unsigned int dim>
void solveSteadyStateProblem(const Config& config);
template<unsigned int dim>
void solveTimeDependentProblem(const Config& config);

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.toml>" << std::endl;
        return 1;
    }

#ifdef _OPENMP
    //Output OpenMP info (max threads, parameters, etc.)
    int max_threads = omp_get_max_threads();
    std::cout << "OpenMP is enabled." << std::endl;
    std::cout << "Max threads: " << max_threads << std::endl;
#else
    std::cout << "OpenMP is not enabled. Running sequentially." << std::endl;
#endif

    try {
        Config config = Config::loadFromFile(argv[1]);

        // Check if problem is time-dependent
        if (config.problem.time_dependent) {
            std::cout << "=== TIME-DEPENDENT PROBLEM ===" << std::endl;
            switch(config.problem.dimension) {
                case 1:
                    solveTimeDependentProblem<1>(config);
                    break;
                case 2:
                    solveTimeDependentProblem<2>(config);
                    break;
                case 3:
                    solveTimeDependentProblem<3>(config);
                    break;
                default:
                    std::cerr << "Error: Unsupported dimension " << config.problem.dimension << std::endl;
                    exit(-1);
            }
        } else {
            std::cout << "=== STEADY-STATE PROBLEM ===" << std::endl;
            switch(config.problem.dimension) {
                case 1:
                    solveSteadyStateProblem<1>(config);
                    break;
                case 2:
                    solveSteadyStateProblem<2>(config);
                    break;
                case 3:
                    solveSteadyStateProblem<3>(config);
                    break;
                default:
                    std::cerr << "Error: Unsupported dimension " << config.problem.dimension << std::endl;
                    exit(-1);
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

template<unsigned int dim>
void solveSteadyStateProblem(const Config& config) {
    if (dim < 1 || dim > 3) {
        std::cerr << "Error: Unsupported dimension " << dim << std::endl;
        exit(-1);
    } 

    Grid<dim> grid = config.createGrid<dim>();
    Function<dim,1> forcing = config.createForcingFunction<dim>();
    Function<dim,1> diffusion = config.createDiffusionFunction<dim>();
    Function<dim,dim> transport = config.createTransportFunction<dim>();
    Function<dim,1> reaction = config.createReactionFunction<dim>();
    BoundaryConditions<dim,1> bc = config.createBoundaryConditions<dim>();
    std::unique_ptr<QuadratureRule<dim>> quadrature = config.createQuadrature<dim>();

    Fem<dim> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
    fem.assemble();
    fem.solve();
    fem.outputCsv(config.problem.output_file + ".csv");
    fem.outputVtu(config.problem.output_file + ".vtu");
}

template<unsigned int dim>
void solveTimeDependentProblem(const Config& config) {
    if (dim < 1 || dim > 3) {
        std::cerr << "Error: Unsupported dimension " << dim << std::endl;
        exit(-1);
    }

    Grid<dim> grid = config.createGrid<dim>();
    Function<dim,1> diffusion = config.createDiffusionFunction<dim>();
    Function<dim,dim> transport = config.createTransportFunction<dim>();
    Function<dim,1> reaction = config.createReactionFunction<dim>();
    BoundaryConditions_td<dim,1> bc = config.createBoundaryConditionsTD<dim>();
    std::unique_ptr<QuadratureRule<dim>> quadrature = config.createQuadrature<dim>();
    
    FemTD<dim> femtd(grid, diffusion, transport, reaction, bc, *quadrature);
    
    // Set time-dependent forcing (simplified)
    femtd.set_forcing(
        config.createForcingFunction_td<dim>()
    );
    
    // Set initial condition
    Function<dim,1> initial = config.createInitialConditionFunction<dim>();
    femtd.set_initial_condition(initial);
    
    // Run simulation
    femtd.run(config.time_dependent.final_time, 
                config.time_dependent.time_step, 
                config.time_dependent.theta,
                config.problem.output_file, config.problem.output_file);
                // VTU prefix               // CSV prefix
}
