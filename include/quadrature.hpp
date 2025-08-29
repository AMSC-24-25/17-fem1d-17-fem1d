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
        std::vector<double>& weights
    ) const;
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

// Class for quadrature on boundary segments (Gauss-Legendre)
template<unsigned int dim>
class GaussLegendre {
private:
    std::vector<std::vector<double>> points;        // Quadrature points 
    std::vector<double> w;             // Quadrature weights
    
public:
    // Constructor for 2-point quadrature (order 3)
    GaussLegendre();
    
    // Main method to obtain quadrature data
    void getQuadratureData(const BoundaryCell<dim>& edge,
                          std::vector<Point<dim+1>>& quadrature_points,
                          std::vector<std::vector<double>>& phi,
                          std::vector<double>& weights) const;

    // Integrates shape functions along a segment for Neumann BC
    void integrateShapeFunctions(const BoundaryCell<dim>& edge,
                                const Function<dim+1,1>& neumannFunc,
                                std::vector<double>& contributions) const;
};

#endif  // QUADRATURE