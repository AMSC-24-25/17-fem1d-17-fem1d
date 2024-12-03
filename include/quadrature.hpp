#ifndef QUADRATURE
#define QUADRATURE

#include "function.hpp"

class QuadratureBase{
    public:
    QuadratureBase(Function _f) : function(_f) {}

    virtual ~QuadratureBase(){}

    virtual double integrate(double a, double b) const = 0;

    protected:

    const Function function; 
};

class SimpsonQuadrature : public QuadratureBase{
    public:
    SimpsonQuadrature(Function _f) : QuadratureBase(_f) {}

    double integrate(double a, double b) const override;

};

class TwoPointsQuadrature : public QuadratureBase{
    public:
    TwoPointsQuadrature(Function _f) : QuadratureBase(_f) {}

    double integrate(double a, double b) const override;

};


#endif  // QUADRATURE