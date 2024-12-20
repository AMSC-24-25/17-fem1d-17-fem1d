#ifndef GRID_1D
#define GRID_1D

#include "function.hpp"

#include <iostream>
#include <vector>

using FunctionVector = std::vector<Function>;

class Grid1D{
    private:
    const double start;
    const double end;

    const int N;
    const double h;

    public:
    Grid1D(double start, double end, int N) :
        start(start), end(end), N(N), h((end-start)/(N-1))
        {}

        
    const double getEnd() const {
        return end;
    }

    double operator()(int k) const{
        if(k > N){
            std::cerr << "k > N passed to Grid\n";
            exit(-1);
        }
        return start + k*h;
    }


    inline double getH() const noexcept { return h; }
    inline int getN() const noexcept { return N; }
    
    FunctionVector getPhiFunctions() const;
};


#endif