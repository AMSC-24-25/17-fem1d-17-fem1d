#ifndef PHI_FUNCTION
#define PHI_FUNCTION

#include "function.hpp"
#include "grid1D.hpp"

class PhiFunction : public Function {
    private:
    int i;
    Grid1D grid;
    
    public:

    PhiFunction(int i, Grid1D _grid) :
        Function(0,0), // doesn't waste resources allocating the lambdas since PhiFunction is already defined
        i(i), grid(_grid) 
    {}

    double value(double x) const override;  
    double grad(double x) const override;  
};


#endif