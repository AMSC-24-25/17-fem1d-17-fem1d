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
inline void localMatricesP1(
    const Cell<2>& cell,
    const std::function<double(const Point<2>&)>& diffusion,
    const std::function<double(const Point<2>&)>& react,
    const std::function<Point<2>(const Point<2>&)>& transport,
    const std::function<double(const Point<2>&)>& forcing,
    Matrix3d &diffusionLocal,
    Matrix3d &transportLocal,
    Matrix3d &reactionLocal,
    Vector3d &forcingLocal
) {
    std::array<Point<2>,3> gradL;
    double area = triangleArea(cell, gradL);
    diffusionLocal.setZero(); transportLocal.setZero(); 
    reactionLocal.setZero(); forcingLocal.setZero();

    QuadRuleTri Q = quadTriOrder2();
    for(size_t q=0;q<Q.barycPoints.size();++q){
        double l1=Q.barycPoints[q][0], l2=Q.barycPoints[q][1], l3=Q.barycPoints[q][2];
        Point<2> x(
            l1*cell[0][0] + l2*cell[1][0] + l3*cell[2][0],
            l1*cell[0][1] + l2*cell[1][1] + l3*cell[2][1]
        );
        double weight = Q.w[q]*area;
        double diff_loc = diffusion(x);
        double react_loc = react(x);
        Point<2> b = transport(x);
        for(int i=0;i<3;++i){
            double phi_i = (i == 0) ? l1 :
                           (i == 1) ? l2 :
                           l3;
            forcingLocal[i] += forcing(x)*phi_i*weight;
            for(int j=0;j<3;++j){
                double phi_j = (j == 0) ? l1 : 
                               (j == 1) ? l2 : 
                               l3;
                double gradDot = (gradL[j][0]*gradL[i][0] + gradL[j][1]*gradL[i][1]);
                diffusionLocal(i,j) += diff_loc * gradDot * weight;
                transportLocal(i,j) += (b[0]*gradL[j][0] + b[1]*gradL[j][1]) * phi_i * weight; // check sign via weak form
                reactionLocal(i,j) += react_loc * phi_i * phi_j * weight;
            }
        }
    }
}

#endif  // QUADRATURE