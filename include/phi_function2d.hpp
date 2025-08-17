// PhiFunction2D.hpp
#ifndef PHI_FUNCTION_2D
#define PHI_FUNCTION_2D

#include "function.hpp"
#include "cell.hpp"
#include "point.hpp"

class PhiFunction2D : public Function {
private:
    int nodeIndex;
    Cell<2> triangle;

public:
    // i: indice del nodo locale (0, 1, 2), cell: triangolo di riferimento
    PhiFunction2D(int i, const Cell<2>& cell);
    
    // Valuta la shape function nel punto p (coordinate baricentriche)
    double value2D(const Point<2>& p) const;
    
    // Gradiente della shape function (costante su triangolo)
    Point<2> gradient2D() const;
};

#endif