#include "quadrature.hpp"

template class QuadratureRule<1>;
template class QuadratureRule<2>;
template class QuadratureRule<3>;

template<>
OrderTwoQuadrature<1>::OrderTwoQuadrature() {
    // Two-point Gauss-Legendre quadrature for interval [-1,1]
    barycPoints = {
        {0.5 + 0.5 * (-1.0/std::sqrt(3.0)), 0.5 - 0.5 * (-1.0/std::sqrt(3.0))},
        {0.5 + 0.5 * (1.0/std::sqrt(3.0)), 0.5 - 0.5 * (1.0/std::sqrt(3.0))}
    };
    w = {0.5, 0.5};
}

template<>
OrderFourQuadrature<1>::OrderFourQuadrature() {
    // Four-point Gauss-Legendre quadrature for interval [-1,1]
    double sqrt_30 = std::sqrt(30.0);
    double x1 = -std::sqrt((3.0 + 2.0 * sqrt_30 / 5.0) / 7.0);
    double x2 = -std::sqrt((3.0 - 2.0 * sqrt_30 / 5.0) / 7.0);
    double x3 = -x2;
    double x4 = -x1;

    double w1 = (18.0 - sqrt_30) / 36.0;
    double w2 = (18.0 + sqrt_30) / 36.0;
    double w3 = w2;
    double w4 = w1;

    barycPoints = {
        {0.5 + 0.5 * x1, 0.5 - 0.5 * x1},
        {0.5 + 0.5 * x2, 0.5 - 0.5 * x2},
        {0.5 + 0.5 * x3, 0.5 - 0.5 * x3},
        {0.5 + 0.5 * x4, 0.5 - 0.5 * x4}
    };
    w = {w1, w2, w3, w4};
}

template<>
OrderTwoQuadrature<3>::OrderTwoQuadrature() {
    // Four-point quadrature rule for tetrahedron (degree 2 exact)
    // Barycentric coordinates and weights:
    barycPoints = {
        {13.0/22.0, 3.0/22.0, 3.0/22.0, 3.0/22.0},
        {3.0/22.0, 13.0/22.0, 3.0/22.0, 3.0/22.0},
        {3.0/22.0, 3.0/22.0, 13.0/22.0, 3.0/22.0},
        {3.0/22.0, 3.0/22.0, 3.0/22.0, 13.0/22.0}
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
    // Six-point quadrature rule for triangle (degree 4 exact, symmetric)
    barycPoints = {
        {5.0/11.0, 5.0/11.0, 1.0/11.0},
        {5.0/11.0, 1.0/11.0, 5.0/11.0},
        {1.0/11.0, 5.0/11.0, 5.0/11.0},
        {1.0/11.0, 1.0/11.0, 9.0/11.0},
        {1.0/11.0, 9.0/11.0, 1.0/11.0},
        {9.0/11.0, 1.0/11.0, 1.0/11.0}
    };
    // Weights sum to 1
    w = {
        2.0/9.0, 2.0/9.0, 2.0/9.0,
        1.0/9.0, 1.0/9.0, 1.0/9.0
    };
}

template<unsigned int dim>
void QuadratureRule<dim>::getQuadratureData(
    const Cell<dim>& cell,
    std::vector<Point<dim>>& grad_phi,
    std::vector<Point<dim>>& quadrature_points,
    std::vector<std::vector<double>>& phi,
    std::vector<double>& weights
){
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
        std::array<double, dim> pCoords;
        for (size_t i = 0; i < dim; i++){
           pCoords[i] = 0;
           for (size_t j = 0; j < dim + 1; j++)
               pCoords[i] += barycPoints[q][j] * cell[j][i];
        }

        quadrature_points.push_back(Point<dim>(pCoords));
        phi.push_back(barycPoints[q]);
        weights.push_back(w[q] * cellMeasure);
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

//TODO non so se csia giusto Function <2,1>
double GaussLegendre1D::integrate(const BoundaryCell<1>& edge, 
                                 const Function<2,1>& func) const {
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


//TODO non so se csia giusto Function <2,1>

void GaussLegendre1D::integrateShapeFunctions(const BoundaryCell<1>& edge,
                                             const Function<2, 1>& neumannFunc,
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

// ============= IMPLEMENTAZIONE GaussLegendre2D =============

void GaussLegendre2D::getQuadratureData(const BoundaryCell<2>& face,
                                        std::vector<Point<3>>& quadrature_points,
                                        std::vector<std::vector<double>>& phi,
                                        std::vector<double>& weights) const {
    // Area della faccia triangolare
    double area = faceArea(face);
    
    // Inizializza i vettori di output
    quadrature_points.clear();
    phi.clear();
    weights.clear();
    quadrature_points.reserve(points.size());
    phi.reserve(points.size());
    weights.reserve(points.size());
    
    // Itera sui punti di quadratura
    for (size_t q = 0; q < points.size(); ++q) {
        // Punto di quadratura globale su faccia
        Point<3> globalPoint = mapToGlobalFace(face, points[q][0], points[q][1]);
        quadrature_points.push_back(globalPoint);
        
        // Shape functions nel punto di quadratura
        std::vector<double> phi_q;
        getShapeFunctions2D(points[q][0], points[q][1], phi_q);
        phi.push_back(phi_q);
        
        // Peso fisico (peso * area della faccia)
        weights.push_back(w[q] * area);
    }
}

double GaussLegendre2D::integrate(const BoundaryCell<2>& face, 
                                 const Function<3,1>& func) const {
    std::vector<Point<3>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    getQuadratureData(face, quadrature_points, phi, weights);
    
    double integral = 0.0;
    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        double funcValue = func.value(quadrature_points[q]);
        integral += weights[q] * funcValue;
    }
    
    return integral;
}

void GaussLegendre2D::integrateShapeFunctions(const BoundaryCell<2>& face,
                                             const Function<3,1>& neumannFunc,
                                             std::vector<double>& contributions) const {
    std::vector<Point<3>> quadrature_points;
    std::vector<std::vector<double>> phi;
    std::vector<double> weights;
    
    getQuadratureData(face, quadrature_points, phi, weights);
    
    contributions.resize(3, 0.0);  // 3 nodi per faccia triangolare
    
    for (size_t q = 0; q < quadrature_points.size(); ++q) {
        // Valore della funzione di Neumann nel punto
        double neumannValue = neumannFunc.value(quadrature_points[q]);
        
        // Contributi ai nodi della faccia
        for (int i = 0; i < 3; ++i) {
            contributions[i] += weights[q] * neumannValue * phi[q][i];
        }
    }
}
