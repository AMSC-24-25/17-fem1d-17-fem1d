#include <iostream>
#include "fem1d.hpp"
#include "../include/quadrature.hpp"
double identity(double x){
    return x;
}
 
int main(int argc, char *argv[])
{
    std::cout << "-------------17-FEM1D PROJECT-----------" << std::endl;
    Function function(identity, identity);
    SimpsonQuadrature quadrature(function);
    std::cout << quadrature.integrate(0.0, 10.0) << std::endl;

    return 0;
}