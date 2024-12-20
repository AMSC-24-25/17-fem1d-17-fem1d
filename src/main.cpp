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

    Function forcing(
        [](double x) -> double { //value
            // return sin(2*PI*x);
            return -1;
        }
    );

    Function diffusion_term = OneFunction();
    Function transport_term = Function(
        [](double x) -> double {
            return 0;
        }
    );
    Function reaction_term = Function(
        [](double x) -> double {
            return 0;
        }
    );
    
    bool isNeumann1 = false;
    Function boundary1 = ZeroFunction();
    bool isNeumann2 = true;
    Function boundary2 = ZeroFunction();

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