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

    return ((b-a)/2)*(function(((a+b)/2)) + ((b-a)/2)*(-(1/sqrt(3))) + function(((a+b)/2)) + ((b-a)/2)*(1/sqrt(3)));
}

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
){
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
