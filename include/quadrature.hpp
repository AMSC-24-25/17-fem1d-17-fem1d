#ifndef QUADRATURE
#define QUADRATURE

#include <math.h>
#include "function.hpp"

//da togliere prima o poi
const int dim = 1;

class QuadratureBase{
    public:
    QuadratureBase(Function<dim> f) : function(f) {}

    virtual ~QuadratureBase() = default;

    virtual double integrate(double a, double b) const = 0;

    protected:

    const Function<dim> function; 
};

class MidPointQuadrature : public QuadratureBase{
    public:
    MidPointQuadrature(Function<dim> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

};

class TrapezoidalQuadrature : public QuadratureBase{
    public:
    TrapezoidalQuadrature(Function<dim> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

};

class SimpsonQuadrature : public QuadratureBase{
    public:
    SimpsonQuadrature(Function<dim> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

};

class TwoPointsQuadrature : public QuadratureBase{
    public:
    TwoPointsQuadrature(Function<dim> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

    virtual ~TwoPointsQuadrature() = default;

};


#endif  // QUADRATURE