#include "quadrature.hpp"

void orderTwoQuadrature::getQuadratureData(
    const Cell<2>& cell,
    std::vector<Point<2>>& grad_phi,
    std::vector<Point<2>>& quadrature_points,
    std::vector<std::vector<double>>& phi,
    std::vector<double>& weights


){
    // compute constant gradients of P1 basis functions on this cell
    grad_phi.resize(3);
    double area = cellArea(cell);
    for(int i=0;i<3;++i) grad_phi[i] = barycentricGradient(cell,i);



    // Compute the local matrices by iterating through the quadrature points
    for(size_t q=0; q < barycPoints.size(); ++q){
        Point<3> barycPoint = Point<3>(barycPoints[q]);
        Point<2> p(
            barycPoint[0]*cell[0][0] + barycPoint[1]*cell[1][0] + barycPoint[2]*cell[2][0],
            barycPoint[0]*cell[0][1] + barycPoint[1]*cell[1][1] + barycPoint[2]*cell[2][1]
        );
        quadrature_points.push_back(p);
        std::vector<double> phi_i;
        phi_i.reserve(3);
        for(int i=0;i<3;++i) phi_i.push_back(barycPoint[i]);
        phi.push_back(phi_i);
        weights.push_back( w[q]*area);
    }
}



void FourPointsQuadrature::getQuadratureData(
    const Cell<2>& cell,
    std::vector<Point<2>>& grad_phi,
    std::vector<Point<2>>& quadrature_points,
    std::vector<std::vector<double>>& phi,
    std::vector<double>& weights
){
    // gradienti P1 costanti e area
    grad_phi.resize(3);
    double area = cellArea(cell);
    for (int i = 0; i < 3; ++i) grad_phi[i] = barycentricGradient(cell, i);

    // punti / pesi fisici e valori φ=λ nei punti
    quadrature_points.clear();
    phi.clear();
    weights.clear();
    quadrature_points.reserve(barycPoints.size());
    phi.reserve(barycPoints.size());
    weights.reserve(barycPoints.size());

    for (size_t q = 0; q < barycPoints.size(); ++q) {
        // λ → x fisico
        Point<3> l = Point<3>(barycPoints[q]);
        Point<2> x(
            l[0]*cell[0][0] + l[1]*cell[1][0] + l[2]*cell[2][0],
            l[0]*cell[0][1] + l[1]*cell[1][1] + l[2]*cell[2][1]
        );
        quadrature_points.push_back(x);

        phi.push_back({ l[0], l[1], l[2] });     // P1: φ_i = λ_i
        weights.push_back(w[q] * area);          // peso fisico
    }
}
