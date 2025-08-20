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

class BarycentricQuadRule {

protected:
    std::vector<std::vector<double>> barycPoints; // barycentric pts
    std::vector<double> w;                 // weights summing to 1
    BarycentricQuadRule() = default;              // impedisce istanziazione diretta


public:
    virtual ~BarycentricQuadRule() = default;

    // Assemble local matrices for variable coefficients
    virtual void getQuadratureData(
        const Cell<2>& cell,
        std::vector<Point<2>>& grad_phi,
        std::vector<Point<2>>& quadrature_points,
        std::vector<std::vector<double>>& phi,
        std::vector<double>& weight
    ) = 0;


};


class orderTwoQuadrature : public BarycentricQuadRule {
    public:
    orderTwoQuadrature() {
        barycPoints = {
            {2/3.0, 1/6.0, 1/6.0},
            {1/6.0, 2/3.0, 1/6.0},
            {1/6.0, 1/6.0, 2/3.0}
        };
        w = {1/3.0, 1/3.0, 1/3.0};
    }

    void getQuadratureData(
        const Cell<2>& cell,
        std::vector<Point<2>>& grad_phi,
        std::vector<Point<2>>& quadrature_points,
        std::vector<std::vector<double>>& phi,
        std::vector<double>& weights
    ) override;

};


class FourPointsQuadrature : public BarycentricQuadRule {
public:
    FourPointsQuadrature() {
        // punti in baricentriche (λ1, λ2, λ3)
        barycPoints = {
            {1.0/3.0, 1.0/3.0, 1.0/3.0},   // baricentro
            {1.0/5.0, 1.0/5.0, 3.0/5.0},
            {1.0/5.0, 3.0/5.0, 1.0/5.0},
            {3.0/5.0, 1.0/5.0, 1.0/5.0}
        };
        // pesi normalizzati a sommare 1 (non 1/2): moltiplicherai poi per area
        w = { -27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0 };
    }

    void getQuadratureData(
        const Cell<2>& cell,
        std::vector<Point<2>>& grad_phi,
        std::vector<Point<2>>& quadrature_points,
        std::vector<std::vector<double>>& phi,
        std::vector<double>& weights
    ) override;
};

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