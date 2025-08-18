#ifndef QUADRATURE
#define QUADRATURE

#include "function.hpp"
#include "cell.hpp"
#include "point.hpp"

#include <array>
#include <math.h>
#include <functional>
#include <Eigen/Dense>

using Eigen::Matrix3d; 
using Eigen::Vector3d; 
using Eigen::Vector2d;

class QuadratureBase{
    public:
    QuadratureBase(Function<1> f) : function(f) {}

    virtual ~QuadratureBase() = default;

    virtual double integrate(double a, double b) const = 0;

    protected:

    const Function<1> function; 
};

class MidPointQuadrature : public QuadratureBase{
    public:
    MidPointQuadrature(Function<1> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;
};

class TrapezoidalQuadrature : public QuadratureBase{
    public:
    TrapezoidalQuadrature(Function<1> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;
};

class SimpsonQuadrature : public QuadratureBase{
    public:
    SimpsonQuadrature(Function<1> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;
};

class TwoPointsQuadrature : public QuadratureBase{
    public:
    TwoPointsQuadrature(Function<1> f) : QuadratureBase(f) {}

    double integrate(double a, double b) const override;

    virtual ~TwoPointsQuadrature() = default;

};

struct BarycentricQuadRule {
    std::vector<std::array<double,3>> barycPoints; // barycentric pts
    std::vector<double> w;                 // weights summing to 1

    // Assemble local matrices for variable coefficients
    void localMatricesP1(
        const Cell<2>& cell,
        const Function<2>& diffusion,
        const Function<2>& react,
        const std::function<Point<2>(const Point<2>&)>& transport,
        const Function<2>& forcing,
        Matrix3d &diffusionLocal,
        Matrix3d &transportLocal,
        Matrix3d &reactionLocal,
        Vector3d &forcingLocal
    );
};

inline BarycentricQuadRule quadTriOrder2() {
    BarycentricQuadRule Q;
    Q.barycPoints = {
        {2/3.0, 1/6.0, 1/6.0},
        {1/6.0, 2/3.0, 1/6.0},
        {1/6.0, 1/6.0, 2/3.0}
    };
    Q.w = {1/3.0, 1/3.0, 1/3.0};
    return Q;
}

// Compute area and gradients
inline double cellArea(const Cell<2>& cell, std::array<Point<2>,3>& gradLambda) {
    for(int i=0;i<3;++i) gradLambda[i] = barycentricGradient(cell,i);
    const Point<2>& A = cell[0]; 
    const Point<2>& B = cell[1]; 
    const Point<2>& C = cell[2];

    double detT = (B[0]-A[0])*(C[1]-A[1]) - (C[0]-A[0])*(B[1]-A[1]);
    return 0.5 * std::abs(detT);
}

#endif  // QUADRATURE