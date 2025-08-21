#include <iostream>
#include <fstream>
#include "fem1d.hpp"
#include "fem.hpp"
#include "grid.hpp"
#include "boundary_conditions.hpp"

using std::cout;
using std::endl;

constexpr double PI = EIGEN_PI;

/**
 * Simple FEM solver for 1D and 2D problems
 * Usage: ./fem 1d L N  or  ./fem 2d mesh.msh
 */
int main(int argc, char *argv[])
{
    cout << "-------------17-FEM1D PROJECT-----------" << endl;
    if(argc < 3){
        cout << "Usage: " << argv[0] << " [1d L N] or [2d mesh.msh] or [3d mesh.msh]" << endl;
        cout << "Examples:" << endl;
        cout << "  " << argv[0] << " 1d 5 20" << endl;
        cout << "  " << argv[0] << " 2d square.msh" << endl;
        cout << "  " << argv[0] << " 3d cube.msh" << endl;
        return -1;
    }

    if (argv[1][0] == '1') {
        // 1D case: fix argument parsing
        if (argc < 4) {
            cout << "1D usage: " << argv[0] << " 1d L N" << endl;
            return -1;
        }
        
        Grid1D grid(0, atof(argv[2]), atoi(argv[3]));

        Function<1,1> forcing(
            [](Point<1> p) -> double { //value
                // return sin(2*PI*x);
                return -1;
            }
        );

        Function<1,1> diffusion_term = OneFunction<1,1>();
        Function<1,1> transport_term = Function<1,1>([](Point<1> p) -> double { return 0; });
        Function<1,1> reaction_term = Function<1,1>([](Point<1> p) -> double {return 0; });
    

        BoundaryConditions<1,1> boundary_conditions;
        boundary_conditions.addDirichlet(0, Point<1>(0.0));
        boundary_conditions.addDirichlet(1, Point<1>(0.0));

        OrderTwoQuadrature<1> quadrature;
        Fem<1> fem(grid, forcing, diffusion_term, transport_term, reaction_term, boundary_conditions, quadrature);

        // const char* solPath = "../sol.csv";
        // std::ofstream fsol(solPath);
        fem.assemble();
        fem.solve();
        //fsol.close();

        std::string csvFilePath = "../sol1d.csv";
        std::string vtuFilePath = "../sol1d.vtu";
        fem.outputCsv(csvFilePath);
        fem.outputVtu(vtuFilePath);

        // cout << "Solution:\n" << fem.getSolution() << endl;

        // system("python ../scripts/plot_sol.py");
    }
    else if (argv[1][0] == '2') {
        // 2D case
        Function<2,1> forcing([](Point<2> p) { 
            return 2.0*(p[0] + p[1]) - 2.0*(p[0]*p[0] + p[1]*p[1]);
        });
        

        BoundaryConditions<2,1> boundary_conditions;
        Function<2,1> diffusion([](Point<2> p) { return 1.0; });
        Function<2,1> reaction([](Point<2> p) { return 0.0; });
        Function<2,2> transport([](Point<2> p) { return Point<2>(0.0, 0.0); });

        // 2. Configurazione delle condizioni al contorno PRIMA del parsing

        OrderTwoQuadrature<2> quadrature;

        // Test Neumann: flusso normale specificato sul lato superiore
        boundary_conditions.addNeumann(3, Function<2,1>([](Point<2> p) { 
            return sin(EIGEN_PI * p[0]);  // Flusso sinusoidale lungo x
        }));
        // Configurazione con mix di Dirichlet e Neumann
        boundary_conditions.addDirichlet(0, Point<1>(0.0));
        boundary_conditions.addDirichlet(1, Point<1>(0.0));
        boundary_conditions.addDirichlet(2, Point<1>(0.0));
        boundary_conditions.addDirichlet(3, Point<1>(0.0));

        // Test delle funzioni in alcuni punti
        Point<2> test_points[3] = {Point<2>(0.5, 0.5), Point<2>(0.2, 0.8), Point<2>(0.7, 0.3)};
        cout << "=== MAIN.CPP FUNCTION VALUES ===" << endl;
        for(int i = 0; i < 3; i++) {
            Point<2> p = test_points[i];
            cout << "Point (" << p[0] << ", " << p[1] << "):" << endl;
            cout << "  forcing = " << forcing.value(p) << endl;
            cout << "  diffusion = " << diffusion.value(p) << endl;
            cout << "  reaction = " << reaction.value(p) << endl;
            cout << "  transport = (" << transport.value(p)[0] << ", " << transport.value(p)[1] << ")" << endl;
        }
        cout << "Boundary conditions:" << endl;
        cout << "  Tag 1: Dirichlet u = 0.0" << endl;
        cout << "  Tag 2: Dirichlet u = 0.0" << endl;
        cout << "  Tag 3: Dirichlet u = 0.0" << endl;
        cout << "  Tag 4: Dirichlet u = 0.0" << endl;

        Grid<2> grid;
        grid.parseFromMsh(argv[2]);
        
        // 4. Crea e risolvi il problema FEM con BoundaryConditions
        Fem<2> fem(grid, forcing, diffusion, transport, reaction, boundary_conditions, quadrature);

        fem.assemble();
        fem.solve();

        std::string csvFilePath = "../sol2d.csv";
        std::string vtuFilePath = "../sol2d.vtu";
        fem.outputCsv(csvFilePath);
        fem.outputVtu(vtuFilePath);
    }
    else if (argv[1][0] == '3'){
        // 3D case - Definizione del problema PRIMA del parsing
        cout << "=== Configurazione problema 3D ===" << endl;
        
        // 1. Definizione delle funzioni del problema
        Function<3,1> forcing([](Point<3> p) { 
            return 2 * (p[0] + p[1] + p[2]) - 2 * (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        });
        Function<3,1> diffusion([](Point<3> p) { return 1.0; });
        Function<3,1> reaction([](Point<3> p) { return 0.0; });
        Function<3,3> transport([](Point<3> p) { return Point<3>(0.0, 0.0, 0.0); });

        // 2. Configurazione delle condizioni al contorno PRIMA del parsing
        BoundaryConditions<3,1> boundary_conditions;
        
        // Configurazione con mix di Dirichlet e Neumann
        boundary_conditions.addDirichlet(0, 0.0);
        boundary_conditions.addDirichlet(1, 0.0);
        boundary_conditions.addDirichlet(2, 0.0);
        boundary_conditions.addDirichlet(3, 0.0);

        cout << "Condizioni al contorno configurate:" << endl;
        cout << "  Physical tag 0 (lato sinistro): Dirichlet u = 0" << endl;
        cout << "  Physical tag 1 (lato destro): Dirichlet u = 0" << endl;
        cout << "  Physical tag 2 (lato inferiore): Dirichlet u = 0" << endl;
        cout << "  Physical tag 3 (lato superiore): Dirichlet u = 0" << endl;

        // 3. Parse della mesh
        Grid<3> grid;
        grid.parseFromMsh(argv[2]);
        
        // 4. Crea e risolvi il problema FEM con BoundaryConditions
        Fem<3> fem3d(grid, forcing, diffusion, transport, reaction, boundary_conditions);

        cout << "=== Risoluzione problema FEM 3D ===" << endl;
        std::string csvFilePath = "../sol3d.csv";
        std::string vtuFilePath = "output/sol3d.vtu";
        fem3d.assemble();
        fem3d.solve();
        fem3d.outputCsv(csvFilePath);
        cout << "3D solution saved to " << csvFilePath << endl;
        fem3d.outputVtu(vtuFilePath);
        cout << "3D solution saved to " << vtuFilePath << endl;
    }
    else {
        cout << "First argument must be '1d' or '2d' or '3d'" << endl;
        return -1;
    }

    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'