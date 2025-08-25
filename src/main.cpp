#include <iostream>
#include <fstream>
#include "fem.hpp"
#include "grid.hpp"
#include "grid1D.hpp"
#include "boundary_conditions.hpp"

using std::cout;
using std::endl;

constexpr double PI = EIGEN_PI;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Simple FEM solver for 1D and 2D problems
 * Usage: ./fem 1d L N  or  ./fem 2d mesh.msh
 */


int main(int argc, char *argv[])
{
    cout << "-------------17-FEM1D PROJECT-----------" << endl;


    if(argc < 3){
        cout << "Usage: " << argv[0] << " [1d L N] or [2d mesh.msh] or [3d mesh.msh]" << endl;
        cout << "Optionally, pass num_threads as last argument." << endl;
        cout << "Examples:" << endl;
        cout << "  " << argv[0] << " 1d 5 20" << endl;
        cout << "  " << argv[0] << " 2d square.msh" << endl;
        cout << "  " << argv[0] << " 3d cube.msh" << endl;
        cout << "  " << argv[0] << " 1d 5 20 4       (4 threads)" << endl;
        cout << "  " << argv[0] << " 3d cube.msh 8   (8 threads)" << endl;
        return -1;
    }

    
    #ifdef _OPENMP
        omp_set_dynamic(0);
        if (argv[1][0] == '1'){
            if (argc >= 5) {
                int nThreads = std::atoi(argv[4]);
                omp_set_num_threads(nThreads);
            }
        } else {
            if (argc >= 4){
                int nThreads = std::atoi(argv[3]);
                omp_set_num_threads(nThreads);
            } 
        }

        std::cout << "OpenMP is enabled." << std::endl;
        std::cout << "Max threads: " << omp_get_max_threads() << std::endl;
    #else
        std::cout << "OpenMP is not enabled. Running sequentially." << std::endl;
    #endif

    if (argv[1][0] == '1') {
    // 1D case: fix argument parsing
        if (argc < 4) {
            cout << "1D usage: " << argv[0] << " 1d L N" << endl;
            return -1;
        }

        Grid1D grid(0, atof(argv[2]), atoi(argv[3]));

        Function<1,1> forcing(
            [](Point<1> p) -> double { //value
                // return sin(2*PI*p[0]);
                return (4.0 * M_PI * M_PI * 2.0 + 5.0) * sin(2.0 * M_PI * p[0]) 
                + 2.0 * M_PI * 4.0 * cos(2.0 * M_PI * p[0]);
            }
        );

    // Function<1,1> diffusion_term = OneFunction<1,1>();
        Function<1,1> diffusion_term = Function<1,1>([](Point<1> p) -> double { return 2.0; });
        Function<1,1> transport_term = Function<1,1>([](Point<1> p) -> double { return 4.0; });
        Function<1,1> reaction_term = Function<1,1>([](Point<1> p) -> double {return 5.0; });
    

        BoundaryConditions<1,1> boundary_conditions;
        boundary_conditions.addDirichlet(0, Point<1>(0.0));
        boundary_conditions.addNeumann(1, Point<1>(2.0 * M_PI * 2.0));

        OrderTwoQuadrature<1> quadrature;
        Fem<1> fem(grid, forcing, diffusion_term, transport_term, reaction_term, boundary_conditions, quadrature);

        fem.assemble();
        fem.solve();

        std::string csvFilePath = "output/sol1d.csv";
        std::string vtuFilePath = "output/sol1d.vtu";
        fem.outputCsv(csvFilePath);
        fem.outputVtu(vtuFilePath);

    // system("python ../scripts/plot_sol.py");
    }
    else if (argv[1][0] == '2') {
        // 2D case

        Function<2,1> forcing([](Point<2> p) { 
            return (-1.0*p[1]) + (-1.0*p[0]) + 1.0*p[0] * p[1];
        });
        
        BoundaryConditions<2,1> boundary_conditions;
        Function<2,1> diffusion([](Point<2> p) { return 1.0; });
        Function<2,1> reaction([](Point<2> p) { return 1.0; });
        Function<2,2> transport([](Point<2> p) { return Point<2>(-1.0, -1.0); });

        OrderTwoQuadrature<2> quadrature;

        // Mix of Dirichlet and Neumann
        boundary_conditions.addDirichlet(0, Function<2,1>([](Point<2> p) { return p[0] * p[1]; }));
        boundary_conditions.addNeumann(1, Function<2,1>([](Point<2> p) { return p[1]; }));
        boundary_conditions.addDirichlet(2, Function<2,1>([](Point<2> p) { return p[0] * p[1]; }));
        boundary_conditions.addNeumann(3, Function<2,1>([](Point<2> p) { return p[0]; }));

        cout << "Boundary conditions:" << endl;
        cout << "  Tag 0: Dirichlet u = 0.0" << endl;
        cout << "  Tag 1: Dirichlet u = 0.0" << endl;
        cout << "  Tag 2: Dirichlet u = 0.0" << endl;
        cout << "  Tag 3: Dirichlet u = 0.0" << endl;

        Grid<2> grid;
        grid.parseFromMsh(argv[2]);
        
        Fem<2> fem(grid, forcing, diffusion, transport, reaction, boundary_conditions, quadrature);

        fem.assemble();
        fem.solve();

        std::string csvFilePath = "output/sol2d.csv";
        std::string vtuFilePath = "output/sol2d.vtu";
        fem.outputCsv(csvFilePath);
        fem.outputVtu(vtuFilePath);
    }
    else if (argv[1][0] == '3'){
        cout << "=== 3D Problem Configuration ===" << endl;
        
        Function<3,1> forcing([](Point<3> p) { 
            return -1.0*p[1]*p[2] -1.0*p[0]*p[2] -1.0*p[0]*p[1] + 1.0*p[0] * p[1] * p[2];
        });
        Function<3,1> diffusion([](Point<3> p) { return 1.0; });
        Function<3,1> reaction([](Point<3> p) { return 1.0; });
        Function<3,3> transport([](Point<3> p) { return Point<3>(-1.0, -1.0, -1.0); });

    // 2. Configure boundary conditions BEFORE mesh parsing
        BoundaryConditions<3,1> boundary_conditions;

        boundary_conditions.addDirichlet(0, Function<3,1>([](Point<3> p) { return p[0]*p[1]*p[2]; }));
        boundary_conditions.addNeumann(1, Function<3,1>([](Point<3> p) { return p[1]*p[2]; }));
        // boundary_conditions.addDirichlet(1, Function<3,1>([](Point<3> p) { return p[0]*p[1]*p[2]; }));
        boundary_conditions.addDirichlet(2, Function<3,1>([](Point<3> p) { return p[0]*p[1]*p[2]; }));
        // boundary_conditions.addNeumann(3, Function<3,1>([](Point<3> p) { return p[0]*p[2]; }));
        boundary_conditions.addDirichlet(3, Function<3,1>([](Point<3> p) { return p[0]*p[1]*p[2]; }));
        boundary_conditions.addDirichlet(4, Function<3,1>([](Point<3> p) { return p[0]*p[1]*p[2]; }));
        boundary_conditions.addDirichlet(5, Function<3,1>([](Point<3> p) { return p[0]*p[1]*p[2]; }));

        cout << "Boundary conditions configured:" << endl;
        cout << "  Physical tag 0 (back face): Dirichlet u = xyz" << endl;
        cout << "  Physical tag 1 (front face): Neumann g = yz" << endl;
        cout << "  Physical tag 2 (left face): Dirichlet u = xyz" << endl;
        cout << "  Physical tag 3 (right face): Neumann g = xz" << endl;
        cout << "  Physical tag 4 (bottom face): Dirichlet u = xyz" << endl;
        cout << "  Physical tag 5 (top face): Dirichlet u = xyz" << endl;
        cout << "NOTE: 'front' is considered as the face whose normal is (1, 0, 0)." << endl;

    // 3. Mesh parsing
        Grid<3> grid;
        grid.parseFromMsh(argv[2]);

        cout << "Mesh parsed successfully:" << endl;
        cout << "  Number of elements: " << grid.getNumElements() << endl;
        cout << "  Number of nodes: " << grid.getNumNodes() << endl;

        Fem<3> fem3d(grid, forcing, diffusion, transport, reaction, boundary_conditions, OrderTwoQuadrature<3>());

        cout << "=== Solving 3D FEM problem ===" << endl;
        std::string csvFilePath = "output/sol3d.csv";
        std::string vtuFilePath = "output/sol3d.vtu";
        fem3d.assemble();
        fem3d.solve();
        fem3d.outputCsv(csvFilePath);

        cout << "3D solution saved to " << csvFilePath << endl;
        fem3d.outputVtu(vtuFilePath);
        cout << "3D solution saved to " << vtuFilePath << endl;
    } else {
        cout << "First argument must be '1d' or '2d' or '3d'" << endl;
        return -1;
    }

    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'