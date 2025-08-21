#include "config.hpp"
#include "fem.hpp"
#include <iostream>
#include <filesystem>
#include <fstream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.toml>" << std::endl;
        return 1;
    }
    
    try {
        Config config = Config::loadFromFile(argv[1]);
        std::filesystem::create_directories("../output");
        
        if (config.problem.dimension == 1) {
            Grid1D grid = config.createGrid1D();
            auto forcing = config.createForcingFunction<1>();
            auto diffusion = config.createDiffusionFunction<1>();
            auto transport = config.createTransportFunction1D();
            auto reaction = config.createReactionFunction<1>();
            auto bc = config.createBoundaryConditions<1>();
            auto quadrature = config.createQuadrature<1>();

            Fem<1> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
            fem.solve();
            
            std::string filename = "../" + config.problem.output_file + ".csv";
            std::ofstream file(filename);
            if (file.is_open()) {
                file << "x,solution\n";
                auto solution = fem.getSolution();
                for (int i = 0; i < grid.getN(); ++i) {
                    double x = grid(i);
                    file << x << "," << solution[i] << "\n";
                }
                file.close();
            }
            
        } else if (config.problem.dimension == 2) {
            Grid2D grid = config.createGrid2D();
            auto forcing = config.createForcingFunction<2>();
            auto diffusion = config.createDiffusionFunction<2>();
            auto transport = config.createTransportFunction2D();
            auto reaction = config.createReactionFunction<2>();
            auto bc = config.createBoundaryConditions<2>();
            auto quadrature = config.createQuadrature<2>();

            // Test delle funzioni in alcuni punti
            Point<2> test_points[3] = {Point<2>(0.5, 0.5), Point<2>(0.2, 0.8), Point<2>(0.7, 0.3)};
            std::cout << "=== TOML FUNCTION VALUES ===" << std::endl;
            for(int i = 0; i < 3; i++) {
                Point<2> p = test_points[i];
                std::cout << "Point (" << p[0] << ", " << p[1] << "):" << std::endl;
                std::cout << "  forcing = " << forcing.value(p) << std::endl;
                std::cout << "  diffusion = " << diffusion.value(p) << std::endl;
                std::cout << "  reaction = " << reaction.value(p) << std::endl;
                std::cout << "  transport = (" << transport.value(p)[0] << ", " << transport.value(p)[1] << ")" << std::endl;
            }
            
            // Print boundary conditions
            std::cout << "Boundary conditions:" << std::endl;
            for (const auto& boundary : config.boundary_conditions) {
                std::cout << "  Tag " << boundary.tag << ": " << boundary.type << " = " << boundary.function << std::endl;
            }


            Fem<2> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
            fem.assemble();
            fem.solve();
            fem.outputVtu("../" + config.problem.output_file + ".vtu");
            fem.outputCsv("../" + config.problem.output_file + ".csv");
            
        } else {
            std::cerr << "Error: Unsupported dimension " << config.problem.dimension << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
