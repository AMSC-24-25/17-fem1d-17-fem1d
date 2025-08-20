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

// ============= IMPLEMENTAZIONE GaussLegendre1D =============

void GaussLegendre1D::getQuadratureData(const BoundaryCell<1>& edge,
                                        std::vector<Point<2>>& quadrature_points,
                                        std::vector<std::vector<double>>& phi,
                                        std::vector<double>& weights) const {
    // Proprietà geometriche dell'edge
    double length = edge[0].distance(edge[1]);
    double jacobian = length / 2.0;  // |J| = length/2 per trasformazione [-1,1] -> edge
    
    // Inizializza i vettori di output
    quadrature_points.clear();
    phi.clear();
    weights.clear();
    quadrature_points.reserve(points.size());
    phi.reserve(points.size());
    weights.reserve(points.size());
    
    // Itera sui punti di quadratura
    for (size_t q = 0; q < points.size(); ++q) {
        // Punto di quadratura globale
        Point<2> globalPoint = mapToGlobalEdge(edge, points[q]);
        quadrature_points.push_back(globalPoint);
        
        // Shape functions nel punto di quadratura
        std::vector<double> phi_q;
        getShapeFunctions1D(points[q], phi_q);
        phi.push_back(phi_q);
        
        // Peso fisico
        weights.push_back(w[q] * jacobian);
    }
}

double GaussLegendre1D::integrate(const BoundaryCell<1>& edge, 
                                 const Function<2>& func) const {
    std::vector<Point<2>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    getQuadratureData(edge, quadrature_points, phi, weights);
    
    double integral = 0.0;
    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        double funcValue = func.value(quadrature_points[q]);
        integral += weights[q] * funcValue;
    }
    
    return integral;
}

void GaussLegendre1D::integrateShapeFunctions(const BoundaryCell<1>& edge,
                                             const Function<2>& neumannFunc,
                                             std::vector<double>& contributions) const {
    std::vector<Point<2>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    getQuadratureData(edge, quadrature_points, phi, weights);
    
    contributions.resize(2, 0.0);  // 2 nodi per edge
    
    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        // Valore della funzione di Neumann nel punto
        double neumannValue = neumannFunc.value(quadrature_points[q]);
        
        // Contributi ai nodi dell'edge
        for (int i = 0; i < 2; ++i) {
            contributions[i] += weights[q] * neumannValue * phi[q][i];
        }
    }
}
