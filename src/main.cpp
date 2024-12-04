#include <iostream>
#include <math.h>
#include "fem1d.hpp"
 
using std::cout;
using std::endl;

constexpr double PI = 3.14;

/**
 * Arg1 = start
 * Arg2 = end
 * Arg3 = N
 */
int main(int argc, char *argv[])
{
    std::cout << "-------------17-FEM1D PROJECT-----------" << std::endl;
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

    Fem1D fem(
        grid, forcing, ZeroFunction(), diffusion_term, 
        false, false, ZeroFunction(), ZeroFunction()
    );

    fem.assemble();
    fem.solve();

    return 0;
}