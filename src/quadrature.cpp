#include <math.h>
#include "../include/quadrature.hpp"

double MidPointQuadrature::integrate(double a, double b) const{

    return (b-a) * function((a+b)/2);
}

double TrapezoidalQuadrature::integrate(double a, double b) const{

    return ((b-a)/2) * (function(a) + function(b));
}

double SimpsonQuadrature::integrate(double a, double b) const {
    
    return ((b-a)/6) * (function(a) + 4*function((a+b)/2) + function(b));
}

double TwoPointsQuadrature::integrate(double a, double b) const {

    return ((b-a)/2) * (function((a+b)/2) + ((b-a)/2)*(-(1/sqrt(3))) + function(((a+b)/2)) + ((b-a)/2)*(1/sqrt(3)));
}