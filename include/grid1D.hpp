#ifndef GRID_1D
#define GRID_1D

#include <iostream>

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
};


#endif