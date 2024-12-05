#include <iostream>
#include <math.h>
#include "fem1d.hpp"
 
using std::cout;
using std::endl;

constexpr double PI = 3.14;

/**
 * Arg1 = L (end of domain)
 * Arg3 = N (# DoFs)
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
            return (2*PI*x);
        }
    );

    Function diffusion_term = OneFunction();

    Fem1D fem(
        grid, forcing, ZeroFunction(), diffusion_term, 
        false, false, ZeroFunction(), ZeroFunction()
    );

    fem.assemble();
    fem.solve();

    cout << "Solution:\n" << fem.getSolution() << endl;

    return 0;
}

//alias v='valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind-out.txt ./fem1d.exe 5 20'