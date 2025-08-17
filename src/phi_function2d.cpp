// PhiFunction2D.cpp
#include "phi_function2d.hpp"

PhiFunction2D::PhiFunction2D(int i, const Cell<2>& cell)
    : Function(
        // Value: dummy function per compatibilità con Function 1D
        [](Point<2> p) -> double {
            return 0.0; // Non usato in questo contesto 2D
        },
        // Gradient: dummy function per compatibilità
        [](Point<2> p) -> double {
            return 0.0; // Non usato in questo contesto 2D
        }
    ), nodeIndex(i), triangle(cell)
{}

// Valuta la shape function nel punto p usando coordinate baricentriche
double PhiFunction2D::value2D(const Point<2>& p) const {
    auto bary = barycentricCoordinates(triangle, p);
    return bary[nodeIndex];
}

// Restituisce il gradiente della shape function (costante su triangolo)
Point<2> PhiFunction2D::gradient2D() const {
    return barycentricGradient(triangle, nodeIndex);
}