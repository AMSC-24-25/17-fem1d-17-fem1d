#ifndef QUADRATURE
#define QUADRATURE

#include "function.hpp"

class QuadratureBase{
    public:
    QuadratureBase(Function _f) : function(_f) {}

    virtual double integrate(double a, double b) const;

    protected:

    const Function function; 
};

class TrapezoidalQuadrature : public QuadratureBase{
    public:
    TrapezoidalQuadrature(Function _f) : QuadratureBase(_f) {}

    double integrate(double a, double b) const override;

    private: 
    int weightA = 1;
    int weightB = 1;
};


#endif  // QUADRATURE