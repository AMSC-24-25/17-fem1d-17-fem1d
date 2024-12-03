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

class TrapezoidalQuadrature : public QuadratureBase{
    public:
    //virtual ~TrapezoidalQuadrature() = default;
    TrapezoidalQuadrature(Function _f) : QuadratureBase(_f) {}

    double integrate(double a, double b) const override;

};


#endif  // QUADRATURE