#ifndef QUADRATURE
#define QUADRATURE

#include <math.h>
#include "function.hpp"

class QuadratureBase{
    public:
    QuadratureBase(Function f) : function(f) {}

    virtual ~QuadratureBase() = default;

    virtual double integrate(double a, double b) const = 0;

    protected:

    const Function function; 
};

class MidPointQuadrature : public QuadratureBase{
    public:
    MidPointQuadrature(Function f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

};

class TrapezoidalQuadrature : public QuadratureBase{
    public:
    TrapezoidalQuadrature(Function f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

};

class SimpsonQuadrature : public QuadratureBase{
    public:
    SimpsonQuadrature(Function f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

};

class TwoPointsQuadrature : public QuadratureBase{
    public:
    TwoPointsQuadrature(Function f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

    virtual ~TwoPointsQuadrature() = default;

};


#endif  // QUADRATURE