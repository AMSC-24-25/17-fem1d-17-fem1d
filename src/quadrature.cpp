#include "../include/quadrature.hpp"
#include "../include/phi_function2d.hpp"

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

void QuadRuleTri::localMatricesP1(
    const Cell<2>& cell,
    const Function<2>& diffusion,
    const Function<2>& react,
    const std::function<Point<2>(const Point<2>&)>& transport,
    const Function<2>& forcing,
    Matrix3d &diffusionLocal,
    Matrix3d &transportLocal,
    Matrix3d &reactionLocal,
    Vector3d &forcingLocal
){
    std::array<Point<2>,3> grad_phi;
    double area = cellArea(cell, grad_phi);
    diffusionLocal.setZero(); transportLocal.setZero(); 
    reactionLocal.setZero(); forcingLocal.setZero();

    for(size_t q=0; q < barycPoints.size(); ++q){
        Point<3> barycPoint = Point<3>(barycPoints[q]);
        Point<2> p(
            barycPoint[0]*cell[0][0] + barycPoint[1]*cell[1][0] + barycPoint[2]*cell[2][0],
            barycPoint[0]*cell[0][1] + barycPoint[1]*cell[1][1] + barycPoint[2]*cell[2][1]
        );
        double weight = w[q]*area;
        double diff_loc = diffusion(p);
        double react_loc = react(p);
        Point<2> b = transport(p);
        for(int i=0;i<3;++i){
            double phi_i = barycPoint[i];
            // LHS
            for(int j=0;j<3;++j){
                double phi_j = barycPoint[j];
                double gradDot = (grad_phi[j][0]*grad_phi[i][0] + grad_phi[j][1]*grad_phi[i][1]);
                double bDot = (b[0]*grad_phi[j][0] + b[1]*grad_phi[j][1]);
                diffusionLocal(i,j) += diff_loc * gradDot * weight;
                transportLocal(i,j) += bDot * phi_i * weight; // check sign via weak form
                reactionLocal(i,j) += react_loc * phi_i * phi_j * weight;
            }
            
            // RHS
            forcingLocal[i] += forcing(p)*phi_i*weight;
        }
    }
}
