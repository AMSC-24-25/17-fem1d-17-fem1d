/**
 * @file quadrature.hpp
 * @brief Quadrature rules for numerical integration over finite elements (1D, 2D, 3D).
 *
 * This header defines classes for quadrature rules on simplices (segments, triangles, tetrahedra)
 * used in finite element methods. Quadrature rules are provided for:
 *   - 1D segments (Gauss-Legendre)
 *   - 2D triangles (order-2, order-4)
 *   - 3D tetrahedra (order-2, order-4)
 */
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

template<unsigned int dim>
class QuadratureRule {

protected:
    std::vector<std::vector<double>> barycPoints; // barycentric pts
    std::vector<double> w;                 // weights summing to 1
    QuadratureRule() = default;              // prevents direct instantiation

public:
    virtual ~QuadratureRule() = default;

    // Assemble local matrices for variable coefficients
    void getQuadratureData(
        const Cell<dim>& cell,
        std::vector<Point<dim>>& grad_phi,
        std::vector<Point<dim>>& quadrature_points,
        std::vector<std::vector<double>>& phi,
        std::vector<double>& weight
    );
};

template<unsigned int dim>
class OrderTwoQuadrature : public QuadratureRule<dim> {
public:
    OrderTwoQuadrature();
};

template<unsigned int dim>
class OrderFourQuadrature : public QuadratureRule<dim> {
public:
    OrderFourQuadrature();
};

// ============= Boundary Quadrature ===================

// Class for 1D quadrature on boundary segments (Gauss-Legendre)
class GaussLegendre1D {
private:
    std::vector<double> points;        // Quadrature points in [0, 1]
    std::vector<double> w;             // Quadrature weights (same name as other classes)
    
public:
    // Constructor for 2-point quadrature (order 3)
    GaussLegendre1D() {
        // Points and weights for [0, 1]
        double x1 = 0.5 - 0.5/std::sqrt(3.0);
        double x2 = 0.5 + 0.5/std::sqrt(3.0);
        points = {x1, x2};
        w = {0.5, 0.5};
    }
    
    // Main method to obtain quadrature data
    void getQuadratureData(const BoundaryCell<1>& edge,
                          std::vector<Point<2>>& quadrature_points,
                          std::vector<std::vector<double>>& phi,
                          std::vector<double>& weights) const;
    
    // Integrates a function along a boundary segment
    double integrate(const BoundaryCell<1>& edge, 
                    const Function<2,1>& func) const;
    
    // Integrates shape functions along a segment for Neumann BC
    void integrateShapeFunctions(const BoundaryCell<1>& edge,
                                const Function<2,1>& neumannFunc,
                                std::vector<double>& contributions) const;
};

// Class for 2D quadrature on boundary triangles (for 3D)
class GaussLegendre2D {
private:
    std::vector<std::vector<double>> points; // Quadrature points in barycentric coordinates (xi, eta)
    std::vector<double> w;                   // Quadrature weights
    
public:
    // Constructor for 3-point quadrature on triangle (order 2)
    GaussLegendre2D() {
    // Quadrature points on reference triangle (order 2, exact for degree 2 polynomials)
    points = {{1.0/6.0, 1.0/6.0},    // point 1: (1/6, 1/6) - barycentric coordinates (2/3, 1/6, 1/6)
          {2.0/3.0, 1.0/6.0},    // point 2: (2/3, 1/6) - barycentric coordinates (1/6, 2/3, 1/6)
          {1.0/6.0, 2.0/3.0}};   // point 3: (1/6, 2/3) - barycentric coordinates (1/6, 1/6, 2/3)
    w = {1.0/6.0, 1.0/6.0, 1.0/6.0}; // Weights (sum to 1/2, area of reference triangle)
    }
    
    // Main method to obtain quadrature data on a face
    void getQuadratureData(const BoundaryCell<2>& face,
                          std::vector<Point<3>>& quadrature_points,
                          std::vector<std::vector<double>>& phi,
                          std::vector<double>& weights) const;
    
    // Integrates a function on a triangular face
    double integrate(const BoundaryCell<2>& face, 
                    const Function<3,1>& func) const;
    
    // Integrates shape functions on a face for Neumann BC 3D
    void integrateShapeFunctions(const BoundaryCell<2>& face,
                                const Function<3,1>& neumannFunc,
                                std::vector<double>& contributions) const;
};

#endif  // QUADRATURE