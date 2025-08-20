#include <iostream>
#include <fstream>
#include "fem1d.hpp"
#include "fem2d.hpp"
#include "grid2D.hpp"
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
        cout << "Usage: " << argv[0] << " [1d L N] or [2d mesh.msh]" << endl;
        cout << "Examples:" << endl;
        cout << "  " << argv[0] << " 1d 5 20" << endl;
        cout << "  " << argv[0] << " 2d square.msh" << endl;
        return -1;
    }

    if (argv[1][0] == '1') {
        // 1D case: fix argument parsing
        if (argc < 4) {
            cout << "1D usage: " << argv[0] << " 1d L N" << endl;
            return -1;
        }
        
        Grid1D grid(0, atof(argv[2]), atoi(argv[3]));

        Function<1> forcing(
            [](Point<1> p) -> double { //value
                // return sin(2*PI*x);
                return -1;
            }
        );

        Function<1> diffusion_term = OneFunction<1>();
        Function<1> transport_term = Function<1>(
            [](Point<1> p) -> double {
                return 0;
            }
        );
        Function<1> reaction_term = Function<1>(
            [](Point<1> p) -> double {
                return 0;
            }
        );
        
        bool isNeumann1 = false;
        Function<1> boundary1 = ZeroFunction<1>();
        bool isNeumann2 = true;
        Function<1> boundary2 = ZeroFunction<1>();

        Fem1D fem(
            grid, forcing, diffusion_term, transport_term, reaction_term,
            isNeumann1, isNeumann2, boundary1, boundary2
        );
        
        const char* solPath = "../sol.csv";
        std::ofstream fsol(solPath);
        cout << "### ASSEMBLE ###" << endl;
        fem.assemble();
        cout << "### SOLVE    ###" << endl;
        fem.solve(fsol);
        cout << "### DONE     ###" << endl;
        cout << "Solution graph saved in " << solPath << endl;
        fsol.close();

        // cout << "Solution:\n" << fem.getSolution() << endl;

        // system("python ../scripts/plot_sol.py");
    }
    else if (argv[1][0] == '2') {
        // 2D case - Definizione del problema PRIMA del parsing
        cout << "=== Configurazione problema 2D ===" << endl;
        
        // 1. Definizione delle funzioni del problema
        Function<2> forcing([](Point<2> p) { 
            return 2 * (p[0] + p[1]) - 2 * (p[0] * p[0] + p[1] * p[1]);
        });
        Function<2> diffusion([](Point<2> p) { return 1.0; });
        Function<2> reaction([](Point<2> p) { return 0.0; });
        Function<2> transport([](Point<2> p) { return 0.0; });
        
        // 2. Configurazione delle condizioni al contorno PRIMA del parsing
        BoundaryConditions boundaryConditions;
        
        // Configurazione con mix di Dirichlet e Neumann per testare
        boundaryConditions.addDirichlet(0, 0.0);  // Lato sinistro: Dirichlet u = 0
        boundaryConditions.addDirichlet(1, 0.0);  // Lato destro: Dirichlet u = 0  
        boundaryConditions.addDirichlet(2, 0.0);  // Lato inferiore: Dirichlet u = 0
        
        // Test Neumann: flusso normale specificato sul lato superiore
        boundaryConditions.addNeumann(3, Function<2>([](Point<2> p) { 
            return sin(EIGEN_PI * p[0]);  // Flusso sinusoidale lungo x
        }));

        cout << "Condizioni al contorno configurate:" << endl;
        cout << "  Physical tag 0 (lato sinistro): Dirichlet u = 0" << endl;
        cout << "  Physical tag 1 (lato destro): Dirichlet u = 0" << endl;
        cout << "  Physical tag 2 (lato inferiore): Dirichlet u = 0" << endl;
        cout << "  Physical tag 3 (lato superiore): NEUMANN du/dn = sin(π*x)" << endl;

        // 3. Parse della mesh
        Grid2D grid;
        grid.parseFromMsh(argv[2]);
        
        // 4. Crea e risolvi il problema FEM con BoundaryConditions
        Fem2D fem2d(grid, forcing, diffusion, transport, reaction, boundaryConditions);
        
        cout << "=== Risoluzione problema FEM 2D ===" << endl;
        std::ofstream fsol("../sol2d.csv");
        fem2d.assemble();
        fem2d.solve(fsol);
        cout << "2D solution saved to sol2d.csv" << endl;
        fem2d.outputVtk("output/sol2d.vtu");
        cout << "2D solution saved to ./output/output_sol2d.vtu" << endl;
        fsol.close();
    }
    else {
        cout << "First argument must be '1d' or '2d'" << endl;
        return -1;
    }

    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'