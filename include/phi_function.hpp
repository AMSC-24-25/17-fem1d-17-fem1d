#ifndef PHI_FUNCTION
#define PHI_FUNCTION

#include "function.hpp"
#include "grid1D.hpp"

class PhiFunction : public Function<1>{
    public:
    PhiFunction(int i, Grid1D grid);
};


#endif