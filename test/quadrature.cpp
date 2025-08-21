#include "quadrature.hpp"
#include "cell.hpp"
#include "point.hpp"

#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <array>

TEST(QuadratureTest, OrderTwo_WeightsAndLinearExactness) {
    // Triangolo unitario: (0,0), (1,0), (0,1) -> area = 1/2
    Cell<2> tri(std::vector<Point<2>>{ 
        Point<2>(std::vector<double>{0.0, 0.0}), 
        Point<2>(std::vector<double>{1.0, 0.0}), 
        Point<2>(std::vector<double>{0.0, 1.0}) 
    }, std::vector<unsigned int>{0u,1u,2u});

    OrderTwoQuadrature<2> quad;
    std::vector<Point<2>> grad_phi;
    std::vector<Point<2>> qp;
    std::vector<std::vector<double>> phi;
    std::vector<double> w;  

    quad.getQuadratureData(tri, grad_phi, qp, phi, w);

    ASSERT_EQ(w.size(), 3u);
    ASSERT_EQ(qp.size(), 3u);

    // La somma dei pesi deve dare l'area
    double sumw = std::accumulate(w.begin(), w.end(), 0.0);
    EXPECT_NEAR(sumw, 0.5, 1e-12);

    // Esattezza su funzione lineare f(x,y)=x+y -> integrale su triangolo unitario = 1/3
    auto f = [](const Point<2>& P) { return P[0] + P[1]; };
    double approx = 0.0;
    for (size_t i = 0; i < qp.size(); ++i) {
        approx += w[i] * f(qp[i]);
    }
    EXPECT_NEAR(approx, 1.0 / 3.0, 1e-12);
}

TEST(QuadratureTest, OrderFour_WeightsAndLinearExactness) {
    // Triangolo unitario: (0,0), (1,0), (0,1) -> area = 1/2
    Cell<2> tri(std::vector<Point<2>>{ 
        Point<2>(std::vector<double>{0.0, 0.0}), 
        Point<2>(std::vector<double>{1.0, 0.0}), 
        Point<2>(std::vector<double>{0.0, 1.0}) 
    }, std::vector<unsigned int>{0u,1u,2u});

    OrderFourQuadrature<2> quad;
    std::vector<Point<2>> grad_phi;
    std::vector<Point<2>> qp;
    std::vector<std::vector<double>> phi;
    std::vector<double> w;

    quad.getQuadratureData(tri, grad_phi, qp, phi, w);

    ASSERT_EQ(w.size(), 4u);
    ASSERT_EQ(qp.size(), 4u);

    // Somma pesi = area
    double sumw = std::accumulate(w.begin(), w.end(), 0.0);
    EXPECT_NEAR(sumw, 0.5, 1e-12);

    // Esatta su f(x,y)=x+y -> integrale su triangolo unitario = 1/3
    auto f = [](const Point<2>& P) { return P[0] + P[1]; };
    double approx = 0.0;
    for (size_t i = 0; i < qp.size(); ++i) {
        approx += w[i] * f(qp[i]);
    }
    EXPECT_NEAR(approx, 1.0 / 3.0, 1e-12);
}



