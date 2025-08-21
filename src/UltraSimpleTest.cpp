#include "config.hpp"
#include "fem.hpp"
#include <iostream>
#include <filesystem>

int main(int argc, char* argv[]) {
    std::cout << "=== ULTRA-SIMPLE FEM Test ===" << std::endl;
    
    try {
        // Hard-coded simple problem (no TOML, no exprtk)
        std::cout << "Creating Grid2D..." << std::endl;
        Grid2D grid;
        grid.parseFromMsh("../mesh/mesh-square-5_gmsh22.msh");
        std::cout << "Grid created. Nodes: " << grid.getNumNodes() 
                  << ", Elements: " << grid.getNumElements() << std::endl;
        
        // Simple lambda functions (replicating working main exactly)
        std::cout << "Creating functions..." << std::endl;
        Function<2,1> forcing([](Point<2> p) { 
            return 2 * (p[0] + p[1]) - 2 * (p[0] * p[0] + p[1] * p[1]);
        });
        Function<2,1> diffusion([](Point<2> p) { return 1.0; });
        Function<2,1> reaction([](Point<2> p) { return 0.0; });
        Function<2,2> transport([](Point<2> p) { return Point<2>(0.0, 0.0); });
        
        // Boundary conditions (replicating working main exactly)
        std::cout << "Creating boundary conditions..." << std::endl;
        BoundaryConditions<2,1> bc;
        
        // Mix di Dirichlet e Neumann come nel main funzionante
        bc.addDirichlet(0, Point<1>(0.0));
        bc.addDirichlet(1, Point<1>(0.0));
        bc.addDirichlet(2, Point<1>(0.0));
        // Tag 3: Neumann invece di Dirichlet!
        bc.addNeumann(3, Function<2,1>([](Point<2> p) { 
            return sin(EIGEN_PI * p[0]);
        }));

        OrderTwoQuadrature<2> quadrature;
        // Create and solve (EXACT sequence from working main)
        std::cout << "Creating Fem<2>..." << std::endl;


        Fem<2> fem(grid, forcing, diffusion, transport, reaction, bc, quadrature);
        
        std::cout << "=== Starting assemble... ===" << std::endl;
        fem.assemble();
        std::cout << "=== Assemble completed ===" << std::endl;
        
        std::cout << "=== Starting solve... ===" << std::endl;
        fem.solve();
        std::cout << "=== SOLVE COMPLETED! ===" << std::endl;
        
        fem.outputVtu("output/ultra_simple.vtu");
        std::cout << "SUCCESS!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
