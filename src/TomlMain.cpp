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

            //from auto to Function<1,1>
            Function<1,1> forcing = config.createForcingFunction<1>();
            Function<1,1> diffusion = config.createDiffusionFunction<1>();
            Function<1,1> transport = config.createTransportFunction1D();
            Function<1,1> reaction = config.createReactionFunction<1>();

            //from auto to BoundaryConditions<1,1>
            BoundaryConditions<1,1> bc = config.createBoundaryConditions<1>();

            //from auto to std::unique_ptr<QuadratureRule<1>>
            std::unique_ptr<QuadratureRule<1>> quadrature = config.createQuadrature<1>();

            Fem<1> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
            fem.assemble();
            fem.solve();
            
            fem.outputCsv("../" + config.problem.output_file + ".csv");
            fem.outputVtu("../" + config.problem.output_file + ".vtu");

        } else if (config.problem.dimension == 2) {
            Grid2D grid = config.createGrid2D();

            //from auto to Function<2,1>
            Function<2,1> forcing = config.createForcingFunction<2>();
            Function<2,1> diffusion = config.createDiffusionFunction<2>();
            Function<2,2> transport = config.createTransportFunction2D();
            Function<2,1> reaction = config.createReactionFunction<2>();

            //from auto to BoundaryConditions<2,1>
            BoundaryConditions<2,1> bc = config.createBoundaryConditions<2>();

            //from auto to std::unique_ptr<QuadratureRule<2>>
            std::unique_ptr<QuadratureRule<2>> quadrature = config.createQuadrature<2>();

            // Test functions at some points
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
            //from const auto& to const BCConfig&
            for (const BCConfig& boundary : config.boundary_conditions) {
                std::cout << "  Tag " << boundary.tag << ": " << boundary.type << " = " << boundary.function << std::endl;
            }


            Fem<2> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
            fem.assemble();
            fem.solve();
            fem.outputVtu("../" + config.problem.output_file + ".vtu");
            fem.outputCsv("../" + config.problem.output_file + ".csv");
            
        } else if (config.problem.dimension == 3) {
            Grid3D grid = config.createGrid3D();

            //from auto to Function<3,1>
            Function<3,1> forcing = config.createForcingFunction<3>();
            Function<3,1> diffusion = config.createDiffusionFunction<3>();
            Function<3,3> transport = config.createTransportFunction3D();
            Function<3,1> reaction = config.createReactionFunction<3>();

            //from auto to BoundaryConditions<3,1>
            BoundaryConditions<3,1> bc = config.createBoundaryConditions<3>();

            //from auto to std::unique_ptr<QuadratureRule<3>>
            std::unique_ptr<QuadratureRule<3>> quadrature = config.createQuadrature<3>();

            // Test functions at some points
            Point<3> test_points[3] = {Point<3>(0.5, 0.5, 0.5), Point<3>(0.2, 0.8, 0.3), Point<3>(0.7, 0.3, 0.9)};
            std::cout << "=== TOML FUNCTION VALUES ===" << std::endl;
            for(int i = 0; i < 3; i++) {
                Point<3> p = test_points[i];
                std::cout << "Point (" << p[0] << ", " << p[1] << ", " << p[2] << "):" << std::endl;
                std::cout << "  forcing = " << forcing.value(p) << std::endl;
                std::cout << "  diffusion = " << diffusion.value(p) << std::endl;
                std::cout << "  reaction = " << reaction.value(p) << std::endl;
                std::cout << "  transport = (" << transport.value(p)[0] << ", " << transport.value(p)[1] << ", " << transport.value(p)[2] << ")" << std::endl;
            }

            // Print boundary conditions
            //from const auto& to const BCConfig&
            for (const BCConfig& boundary : config.boundary_conditions) {
                std::cout << "  Tag " << boundary.tag << ": " << boundary.type << " = " << boundary.function << std::endl;
            }

            Fem<3> fem(grid, forcing, diffusion, transport, reaction, bc, *quadrature);
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
