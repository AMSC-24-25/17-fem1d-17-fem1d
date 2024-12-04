#include <math.h>
#include "../include/quadrature.hpp"
#include "../include/grid1D.hpp"

double MidPointQuadrature::integrate(double a, double b) const{

    return (b-a) * function((a+b)/2);
}

double TrapezoidalQuadrature::integrate(double a, double b) const{

    return ((b-a)/2) * (function(a) + function(b));
}

double SimpsonQuadrature::integrate(double a, double b) const {
    Grid1D grid1D(a, b, 1000);
    double temp = 0.0;
    for (int i = 0; i <= 1000; ++i) {
        temp += function(grid1D(i-1)) + 4*function((grid1D(i-1) + grid1D(i))/2.0) + function(grid1D(i));
    }

    return (grid1D.getH()/6.0) * temp;
}

double TwoPointsQuadrature::integrate(double a, double b) const {

    return ((b-a)/2)*(function(((a+b)/2)) + ((b-a)/2)*(-(1/sqrt(3))) + function(((a+b)/2)) + ((b-a)/2)*(1/sqrt(3)));
}
