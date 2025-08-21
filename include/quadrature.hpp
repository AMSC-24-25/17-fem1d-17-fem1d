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
    QuadratureRule() = default;              // impedisce istanziazione diretta

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

// Classe per quadratura 1D sui segmenti di bordo (Gauss-Legendre)
class GaussLegendre1D {
private:
    std::vector<double> points;        // Punti di quadratura in [-1, 1]
    std::vector<double> w;             // Pesi di quadratura (stesso nome delle altre classi)
    
public:
    // Costruttore per quadratura a 2 punti (ordine 3)
    GaussLegendre1D() {
        // Punti di Gauss-Legendre a 2 punti in [-1, 1]
        points = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
        w = {1.0, 1.0};
    }
    
    // Metodo principale per ottenere dati di quadratura
    void getQuadratureData(const BoundaryCell<1>& edge,
                          std::vector<Point<2>>& quadrature_points,
                          std::vector<std::vector<double>>& phi,
                          std::vector<double>& weights) const;
    
    // Integra una funzione lungo un segmento di bordo
    double integrate(const BoundaryCell<1>& edge, 
                    const Function<2,1>& func) const;
    
    // Integra le shape functions lungo un segmento per Neumann BC
    void integrateShapeFunctions(const BoundaryCell<1>& edge,
                                const Function<2,1>& neumannFunc,
                                std::vector<double>& contributions) const;
};

#endif  // QUADRATURE