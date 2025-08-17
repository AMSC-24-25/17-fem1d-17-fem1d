#include "phi_function.hpp"

PhiFunction::PhiFunction(int i, Grid1D grid) :
        Function(
            // Value
            [grid, i](Point<1> p) -> double {
                double h = grid.getH();

                if(p[0] <= grid(i-1)) return 0.0;
                if(p[0] < grid(i)) return (p[0] - grid(i-1))/h;
                if(p[0] < grid(i+1)) return 1 - (p[0] - grid(i))/h;
                else return 0.0;
            },

            // Gradient
            [grid, i](Point<1> p) -> double {
                int N = grid.getN();

                if(p[0] <= grid(i-1)) return 0.0;
                if(p[0] < grid(i)) return N;
                if(p[0] < grid(i+1)) return -N;
                else return 0.0;
            }
        )
    {}