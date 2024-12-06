#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "fem1d.hpp"
 
using std::cout;
using std::endl;

constexpr double PI = 3.14;

/**
 * Arg1 = L (end of domain)
 * Arg2 = N (# DoFs)
 */
int main(int argc, char *argv[])
{
    cout << "-------------17-FEM1D PROJECT-----------" << endl;
    if(argc < 3){
        cout << "Usage: ./exe L N" << endl;
        exit(-1);
    }

    Grid1D grid(0, atoi(argv[1]), atoi(argv[2]));

    Function forcing(
        [](double x) -> double { //value
            return sin(2*PI*x);
        }
    );

    Function diffusion_term = OneFunction();
    Function reaction_term = ZeroFunction();
    
    bool isNeumann1 = false;
    Function boundary1 = ZeroFunction();
    bool isNeumann2 = false;
    Function boundary2 = ZeroFunction();

    Fem1D fem(
        grid, forcing, reaction_term, diffusion_term, 
        isNeumann1, isNeumann2, boundary1, boundary2
    );

    std::ofstream fsol("sol.csv");
    fem.assemble();
    fem.solve(fsol);
    fsol.close();

    cout << "Solution:\n" << fem.getSolution() << endl;

    // system("python scripts/plot_sol.py");
    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'