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

#include "cell.hpp"
#include "point.hpp"
#include <array>
#include <functional>
#include <Eigen/Dense>

using Eigen::Matrix3d; 
using Eigen::Vector3d; 
using Eigen::Vector2d;

struct QuadRuleTri {
    std::vector<std::array<double,3>> barycPoints; // barycentric pts
    std::vector<double> w;                 // weights summing to 1
};

inline QuadRuleTri quadTriOrder2() {
    QuadRuleTri Q;
    Q.barycPoints = {
        {2.0/3, 1.0/6, 1.0/6},
        {1.0/6, 2.0/3, 1.0/6},
        {1.0/6, 1.0/6, 2.0/3}
    };
    Q.w = {1.0/3,1.0/3,1.0/3};
    return Q;
}

// Compute area and gradients
inline double triangleArea(const Cell<2>& cell, std::array<Point<2>,3>& gradLambda) {
    for(int i=0;i<3;++i) gradLambda[i] = barycentricGradient(cell,i);
    const auto& A = cell[0]; 
    const auto& B = cell[1]; 
    const auto& C = cell[2];

    double detT = (B[0]-A[0])*(C[1]-A[1]) - (C[0]-A[0])*(B[1]-A[1]);
    return 0.5 * std::abs(detT);
}

// Assemble local matrices for variable coefficients
void localMatricesP1(
    const Cell<2>& cell,
    const std::function<double(const Point<2>&)>& diffusion,
    const std::function<double(const Point<2>&)>& react,
    const std::function<Point<2>(const Point<2>&)>& transport,
    const std::function<double(const Point<2>&)>& forcing,
    Matrix3d &diffusionLocal,
    Matrix3d &transportLocal,
    Matrix3d &reactionLocal,
    Vector3d &forcingLocal
);

#endif  // QUADRATURE