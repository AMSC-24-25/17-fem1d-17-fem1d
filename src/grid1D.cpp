#include "../include/grid1D.hpp"

#include "phi_function.hpp"

FunctionVector Grid1D::getPhiFunctions() const{
    //phi_i = 1 on x_i ; 0 on x_j (j != i)
    FunctionVector phiFunctions;
    
    // Phi0 and PhiN are also defined normally outside of grid to avoid code repetition or exceptions
    // It shouldn't cause problems since the user won't ask values outside of grid
    phiFunctions.reserve(N+1);
    for (int k = 0; k <= N; k++)
    {
        phiFunctions.push_back(PhiFunction(k, *this));
    }
    
    return phiFunctions; 
}