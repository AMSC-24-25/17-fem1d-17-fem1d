#include <iostream>
#include <fstream>
#include "fem1d.hpp"
#include "fem2d.hpp"
#include "grid2D.hpp"

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

        // system("python scripts/plot_sol.py");
    }
    else if (argv[1][0] == '2') {
        // 2D case 
        Grid2D grid;
        grid.parseFromMsh(argv[2]);
        
        Function<2> forcing([](Point<2> p) { 
            return -2.0 * PI * PI * sin(PI * p[0]) * sin(PI * p[1]); 
        });
        Function<2> diffusion([](Point<2> p) { return 1.0; }); // Fix: explicit function instead of OneFunction<2>
        Function<2> reaction([](Point<2> p) { return 0.0; });  // Fix: explicit function instead of ZeroFunction<2>
        Function<2> transport([](Point<2> p) { return 0.0; }); // Fix: scalar transport instead of vector
        Function<2> boundary([](Point<2> p) { return 0.0; });  // Fix: explicit function

        Fem2D fem2d(grid, forcing, diffusion, reaction, transport,
                   false, false, boundary, boundary); // Fix: match constructor signature
        
        std::ofstream fsol("../sol2d.csv");
        fem2d.assemble();
        fem2d.solve(fsol);
        cout << "2D solution saved to sol2d.csv" << endl;
        fsol.close();
    }
    else {
        cout << "First argument must be '1d' or '2d'" << endl;
        return -1;
    }

    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'