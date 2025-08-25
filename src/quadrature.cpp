#include "quadrature.hpp"

template class QuadratureRule<1>;
template class QuadratureRule<2>;
template class QuadratureRule<3>;

template<>
OrderTwoQuadrature<1>::OrderTwoQuadrature() {
    // Two-point Gauss-Legendre quadrature for interval [0,1]
    double x1 = 0.7886751346;
    double x2 = 0.2113248654;
    barycPoints = {
        {1.0 - x1, x1},
        {1.0 - x2, x2}
    };
    w = {0.5, 0.5};
}

template<>
OrderFourQuadrature<1>::OrderFourQuadrature() {
    // Four-point Gauss-Legendre quadrature for interval [0,1]
    double x1 = 0.9305681558;
    double x2 = 0.6699905218;
    double x3 = 0.3300094782;
    double x4 = 0.0694318442;

    double w1 = 0.1739274226;
    double w2 = 0.3260725774;

    barycPoints = {
        {1.0 - x1, x1},
        {1.0 - x2, x2},
        {1.0 - x3, x3},
        {1.0 - x4, x4}
    };
    w = {w1, w2, w2, w1};
}

template<>
OrderTwoQuadrature<3>::OrderTwoQuadrature() {
    // Four-point quadrature rule for tetrahedron (degree 2 exact)
    double sqrt5 = std::sqrt(5.0);
    double a = 0.5854101966; 
    double b = 0.1381966011; 
    
    // Barycentric coordinates and weights:
    barycPoints = {
        {b, a, b, b},
        {a, b, b, b},
        {b, b, b, a},
        {b, b, a, b}
    };
    // Weights sum to 1
    w = {1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0};
}
template<>
OrderTwoQuadrature<2>::OrderTwoQuadrature() {
    // Three-point quadrature rule for triangle (order 2)
    barycPoints = {
        {2/3.0, 1/6.0, 1/6.0},
        {1/6.0, 2/3.0, 1/6.0},
        {1/6.0, 1/6.0, 2/3.0}
    };
    // Weights sum to 1
    w = {1/3.0, 1/3.0, 1/3.0};
}
template<>
OrderFourQuadrature<2>::OrderFourQuadrature() {
    double a = 0.4459484909;
    double b = 0.0915762135;
    double c = 0.1081030182;
    double d = 0.8168475730;
    // Six-point quadrature rule for triangle (degree 4 exact, symmetric)
    barycPoints = {
        {a, a, c}, 
        {b, b, d}, 
        {a, c, a}, 
        {b, d, b}, 
        {c, a, a}, 
        {d, b, b}
    };

    // Weights sum to 1
    double w1 = 0.1116907948;
    double w2 = 0.0549758718;
    w = {   
        w1, w2, w1, 
        w2, w1, w2
    };
    for (size_t i = 0; i < w.size(); ++i) {
        w[i] = w[i] * 2.0; // Area of reference triangle is 1
    }
}

template<unsigned int dim>
void QuadratureRule<dim>::getQuadratureData(
    const Cell<dim>& cell,
    std::vector<Point<dim>>& grad_phi,
    std::vector<Point<dim>>& quadrature_points,
    std::vector<std::vector<double>>& phi,
    std::vector<double>& weights
) const {
    // Compute constant gradients of P1 basis functions on this tetrahedral cell
    grad_phi.resize(dim + 1);
    double cellMeasure = cell.measure();
    for(int i=0; i < dim+1; ++i) grad_phi[i] = cell.barycentricGradient(i);

    // Compute quadrature points, shape functions, and weights
    quadrature_points.clear();
    phi.clear();
    weights.clear();
    quadrature_points.reserve(barycPoints.size());
    phi.reserve(barycPoints.size());
    weights.reserve(barycPoints.size());

    for(size_t q=0; q < barycPoints.size(); ++q){
        Point<dim> globalPoint = cell.mapToGlobal(barycPoints[q]);

        quadrature_points.push_back(globalPoint);
        phi.push_back(barycPoints[q]);
        weights.push_back(w[q] * cellMeasure);
    }
}

// ============= GaussLegendre general implementations ======

template<unsigned int dim>
void GaussLegendre<dim>::getQuadratureData(const BoundaryCell<dim>& boundary,
                                          std::vector<Point<dim+1>>& quadrature_points,
                                          std::vector<std::vector<double>>& phi,
                                          std::vector<double>& weights) const {
    // Geometric properties 
    double measure = boundary.measure();
    
    // Initialize output vectors
    quadrature_points.clear();
    phi.clear();
    weights.clear();
    quadrature_points.reserve(points.size());
    phi.reserve(points.size());
    weights.reserve(points.size());
    
    // Iterate over quadrature points
    for (size_t q = 0; q < points.size(); ++q) {
        Point<dim+1> globalPoint = boundary.mapToGlobal(points[q]);
        quadrature_points.push_back(globalPoint);
        
        // Shape functions at the quadrature point (barycentric coordinates)
        phi.push_back(points[q]);
        
        // Physical weight, scaled by cell measure
        weights.push_back(w[q] * measure);
    }
}

template<unsigned int dim>
double GaussLegendre<dim>::integrate(const BoundaryCell<dim>& edge, 
                                 const Function<dim+1,1>& func) const {
    std::vector<Point<dim+1>> quadrature_points;
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

template<unsigned int dim>
void GaussLegendre<dim>::integrateShapeFunctions(const BoundaryCell<dim>& edge,
                                             const Function<dim+1, 1>& neumannFunc,
                                             std::vector<double>& contributions) const {
    std::vector<Point<dim+1>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    getQuadratureData(edge, quadrature_points, phi, weights);

    contributions.resize(dim+1, 0.0);  // 2 nodes for edge, 3 nodes for face

    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        // Value of the Neumann function at the point
        double neumannValue = neumannFunc.value(quadrature_points[q]);
        
        // Contributions to the edge nodes
        for (int i = 0; i < contributions.size(); ++i) {
            contributions[i] += weights[q] * neumannValue * phi[q][i];
        }
    }
}

// ============= GaussLegendre specialized constructors =============

template<>
GaussLegendre<1>::GaussLegendre() {
    // Four-point Gauss-Legendre quadrature for interval [0,1]
    double x1 = 0.9305681558;
    double x2 = 0.6699905218;
    double x3 = 0.3300094782;
    double x4 = 0.0694318442;

    double w1 = 0.1739274226;
    double w2 = 0.3260725774;

    points = {
        {1.0 - x1, x1},
        {1.0 - x2, x2},
        {1.0 - x3, x3},
        {1.0 - x4, x4}
    };
    w = {w1, w2, w2, w1};
}

template<>
GaussLegendre<2>::GaussLegendre() {
    // Three-point quadrature rule for triangle (order 2)
    points = {
        {2/3.0, 1/6.0, 1/6.0},
        {1/6.0, 2/3.0, 1/6.0},
        {1/6.0, 1/6.0, 2/3.0}
    };
    // Weights sum to 1
    w = {1/3.0, 1/3.0, 1/3.0};
}

template class GaussLegendre<1>;
template class GaussLegendre<2>;
