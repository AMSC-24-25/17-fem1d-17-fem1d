#ifndef PROBLEM_FACTORY_HPP
#define PROBLEM_FACTORY_HPP

#include "config.hpp"
#include "grid1D.hpp"
#include "grid2D.hpp" 
#include "boundary_conditions.hpp"
#include "function.hpp"
#include "fem1d.hpp"
#include "fem2d.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>

// Classe base astratta per problemi di qualsiasi dimensione
class ProblemBase {
public:
    virtual ~ProblemBase() = default;
    virtual void solve() = 0;
    virtual void saveResults() = 0;
    virtual void printInfo() const = 0;
};

// Parser semplificato per funzioni matematiche
template<unsigned int dim>
Function<dim,1> parseSimpleFunction(const std::string& expression) {
    // Parser molto semplificato - supporta solo casi base
    if (expression == "0.0" || expression == "0") {
        return Function<dim,1>([](Point<dim>) { return 0.0; });
    }
    if (expression == "1.0" || expression == "1") {
        return Function<dim,1>([](Point<dim>) { return 1.0; });
    }
    if (expression == "x" && dim >= 1) {
        return Function<dim,1>([](Point<dim> p) { return p[0]; }, 
                               [](Point<dim>) { return 1.0; });
    }
    if (expression == "y" && dim >= 2) {
        return Function<dim,1>([](Point<dim> p) { return p[1]; });
    }
    if (expression == "x*x" && dim >= 1) {
        if constexpr (dim == 1) {
            return Function<dim,1>([](Point<dim> p) { return p[0] * p[0]; },
                                   [](Point<dim> p) { return 2.0 * p[0]; });
        } else {
            return Function<dim,1>([](Point<dim> p) { return p[0] * p[0]; });
        }
    }
    if constexpr (dim >= 2) {
        if (expression == "x*x + y*y") {
            return Function<dim,1>([](Point<dim> p) { return p[0]*p[0] + p[1]*p[1]; },
                                   [](Point<dim> p) { return 2.0 * p[0]; },
                                   [](Point<dim> p) { return 2.0 * p[1]; });
        }
    }
    if (expression == "sin(pi*x)" && dim >= 1) {
        if constexpr (dim == 1) {
            return Function<dim,1>([](Point<dim> p) { return sin(M_PI * p[0]); },
                                   [](Point<dim> p) { return M_PI * cos(M_PI * p[0]); });
        } else {
            return Function<dim,1>([](Point<dim> p) { return sin(M_PI * p[0]); });
        }
    }
    if constexpr (dim >= 2) {
        if (expression == "sin(pi*x)*sin(pi*y)") {
            return Function<dim,1>([](Point<dim> p) { return sin(M_PI * p[0]) * sin(M_PI * p[1]); },
                                   [](Point<dim> p) { return M_PI * cos(M_PI * p[0]) * sin(M_PI * p[1]); },
                                   [](Point<dim> p) { return M_PI * sin(M_PI * p[0]) * cos(M_PI * p[1]); });
        }
    }
    
    // Fallback: try to parse as constant
    try {
        double value = std::stod(expression);
        return Function<dim,1>([value](Point<dim>) { return value; });
    } catch (...) {
        std::cerr << "Warning: Could not parse function '" << expression 
                  << "', using zero function" << std::endl;
        return Function<dim,1>([](Point<dim>) { return 0.0; });
    }
}

// Template specializzato per problemi di dimensione specifica
template<unsigned int dim>
class Problem : public ProblemBase {
private:
    Config config;
    Function<dim,1> forcing;
    std::unique_ptr<Fem1D> fem1d;  // For 1D problems
    std::unique_ptr<Fem<dim>> fem2d;  // For 2D problems
    
public:
    Problem(const Config& cfg) : config(cfg), forcing([](Point<dim>) { return 0.0; }) {
        setupBoundaryConditions();
        setupForcingFunction();
        loadAndTestGrid();
    }
    
private:
    void loadAndTestGrid() {
        if constexpr (dim == 1) {
            // Per 1D creiamo una griglia semplice
            std::cout << "Using simple 1D grid (0 to 1, 10 points)" << std::endl;
        } else if constexpr (dim == 2) {
            // Per 2D carichiamo da file
            Grid2D grid2d;
            grid2d.parseFromMsh(config.problem.mesh_file);
            std::cout << "Loaded 2D grid from: " << config.problem.mesh_file << std::endl;
            std::cout << "Grid has " << grid2d.getNumNodes() << " nodes and " 
                      << grid2d.getNumElements() << " elements" << std::endl;
            
            // Test boundary conditions on real grid
            testBoundaryConditions2D(grid2d);
        } else if constexpr (dim == 3) {
            std::cout << "3D grid loading not yet implemented" << std::endl;
            throw std::runtime_error("3D grid not implemented yet");
        }
    }
    
    void testBoundaryConditions2D(Grid2D& grid) {
        BoundaryConditions<2,1> bc;
        
        // Aggiungi le BC dal config
        for (const auto& bcConfig : config.boundary_conditions) {
            Function<2,1> bcFunction = parseSimpleFunction<2>(bcConfig.function);
            
            if (bcConfig.type == BCConfig::DIRICHLET) {
                bc.addDirichlet(bcConfig.tag, bcFunction);
            } else if (bcConfig.type == BCConfig::NEUMANN) {
                bc.addNeumann(bcConfig.tag, bcFunction);
            }
        }
        
        // Test applicazione BC
        int n = grid.getNumNodes();
        Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
        A.setIdentity();
        Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
        
        bc.apply(grid, A, b);
        std::cout << "Applied boundary conditions to " << n << "x" << n << " system" << std::endl;
    }
    
    void testBoundaryConditions1D() {
        BoundaryConditions<1,1> bc;
        
        // Aggiungi le BC dal config
        for (const auto& bcConfig : config.boundary_conditions) {
            Function<1,1> bcFunction = parseSimpleFunction<1>(bcConfig.function);
            
            if (bcConfig.type == BCConfig::DIRICHLET) {
                bc.addDirichlet(bcConfig.tag, bcFunction);
            } else if (bcConfig.type == BCConfig::NEUMANN) {
                bc.addNeumann(bcConfig.tag, bcFunction);
            }
        }
        
        // Test con griglia 1D semplice
        Grid1D grid(0.0, 1.0, 10);
        int n = grid.getN();
        Eigen::SparseMatrix<double, Eigen::RowMajor> A(n, n);
        A.setIdentity();
        Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
        
        bc.apply(grid, A, b);
        std::cout << "Applied 1D boundary conditions to " << n << "x" << n << " system" << std::endl;
    }
    
    void setupBoundaryConditions() {
        std::cout << "Setting up boundary conditions..." << std::endl;
        for (const auto& bcConfig : config.boundary_conditions) {
            std::cout << "  " << (bcConfig.type == BCConfig::DIRICHLET ? "Dirichlet" : "Neumann")
                      << " BC on tag " << bcConfig.tag << " with function: " << bcConfig.function << std::endl;
        }
    }
    
    void setupForcingFunction() {
        forcing = parseSimpleFunction<dim>(config.equation.forcing_function);
        std::cout << "Set forcing function: " << config.equation.forcing_function << std::endl;
    }
    
    void solve() override {
        std::cout << "Solving " << dim << "D problem..." << std::endl;
        
        std::cout << "Equation coefficients:" << std::endl;
        std::cout << "  Diffusion: " << config.equation.diffusion_coefficient << std::endl;
        std::cout << "  Transport: " << config.equation.transport_coefficient << std::endl;
        std::cout << "  Reaction: " << config.equation.reaction_coefficient << std::endl;
        
        std::cout << "Solver settings:" << std::endl;
        std::cout << "  Method: " << config.solver.method << std::endl;
        std::cout << "  Tolerance: " << config.solver.tolerance << std::endl;
        std::cout << "  Max iterations: " << config.solver.max_iterations << std::endl;
        
        // Test boundary conditions
        if constexpr (dim == 1) {
            testBoundaryConditions1D();
            
            // Create and solve 1D problem
            fem1d = std::make_unique<Fem1D>();
            // Initialize fem1d with grid, boundary conditions, and forcing function
            // fem1d->solve();
        } else if constexpr (dim == 2) {
            // Create and solve 2D problem  
            fem2d = std::make_unique<Fem<dim>>();
            // Initialize fem2d with grid, boundary conditions, and forcing function
            // fem2d->solve();
        }
        
        std::cout << "Solution completed!" << std::endl;
    }
    
    void saveResults() override {
        std::cout << "Saving results to: " << config.problem.output_file << std::endl;
        
        // Create output directory if it doesn't exist
        std::filesystem::create_directories("output");
        
        if constexpr (dim == 1) {
            if (fem1d) {
                std::string csv_filename = config.problem.output_file + ".csv";
                fem1d->outputCsv(csv_filename);
                std::cout << "1D results saved to: " << csv_filename << std::endl;
            } else {
                std::cerr << "Warning: No 1D FEM solver created" << std::endl;
            }
        } else if constexpr (dim == 2) {
            if (fem2d) {
                std::string vtu_filename = config.problem.output_file + ".vtu";
                fem2d->outputVtu(vtu_filename);
                std::cout << "2D results saved to: " << vtu_filename << std::endl;
                
                // Save also in CSV for analysis
                std::string csv_filename = config.problem.output_file + ".csv";
                fem2d->outputCsv(csv_filename);
                std::cout << "2D CSV results saved to: " << csv_filename << std::endl;
            } else {
                std::cerr << "Warning: No 2D FEM solver created" << std::endl;
            }
        }
        
        std::cout << "Results saved!" << std::endl;
    }
    
    void printInfo() const override {
        std::cout << "Problem dimension: " << dim << "D" << std::endl;
        std::cout << "Forcing function: " << config.equation.forcing_function << std::endl;
        std::cout << "Boundary conditions: " << config.boundary_conditions.size() << std::endl;
    }
};

// Factory per creare problemi con dimensione runtime
class ProblemFactory {
public:
    static std::unique_ptr<ProblemBase> create(const Config& config) {
        switch (config.problem.dimension) {
            case 1:
                return std::make_unique<Problem<1>>(config);
            case 2:
                return std::make_unique<Problem<2>>(config);
            case 3:
                return std::make_unique<Problem<3>>(config);
            default:
                throw std::runtime_error("Unsupported dimension: " + 
                                        std::to_string(config.problem.dimension));
        }
    }
};

#endif // PROBLEM_FACTORY_HPP
