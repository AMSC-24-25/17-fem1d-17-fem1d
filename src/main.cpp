#include <iostream>
#include <fstream>
#include "fem1d.hpp"

using std::cout;
using std::endl;

constexpr double PI = EIGEN_PI;

/**
 * Arg1 = L (end of domain)
 * Arg2 = N (# DoFs)
 */
int main(int argc, char *argv[])
{
    cout << "-------------17-FEM1D PROJECT-----------" << endl;
    if(argc < 3){
        cout << "Usage: executable L N" << endl;
        exit(-1);
    }

    Grid1D grid(0, atoi(argv[1]), atoi(argv[2]));

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
    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'