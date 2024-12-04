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
        start(start), end(end), N(N), h((end-start)/N) 
        {}

    double operator()(int k) const{
        if(k > N){
            std::cerr << "k > N passed to Grid\n";
            exit(-1);
        }
        return start + k*h;
    }

<<<<<<< HEAD
    double getH() const { return h;}
    double getN() const { return N;}

=======
    inline double getH() const noexcept { return h; }
    inline int getN() const noexcept { return N; }
>>>>>>> a50ccdf96198f3e717fe1e7b1721b5ae8e2e43e7

    FunctionVector getPhiFunctions() const;
};


#endif