#include "phi_function.hpp"

PhiFunction::PhiFunction(int i, Grid1D grid) :
        Function(
            // Value
            [grid, i](double x) -> double {
                double h = grid.getH();

                if(x <= grid(i-1)) return 0.0;
                if(x < grid(i)) return (x - grid(i-1))/h;
                if(x < grid(i+1)) return 1 - (x - grid(i))/h;
                else return 0.0;
            },

            // Gradient
            [grid, i](double x) -> double {
                int N = grid.getN();

                if(x <= grid(i-1)) return 0.0;
                if(x < grid(i)) return N;
                if(x < grid(i+1)) return -N;
                else return 0.0;
            }
        )
    {}