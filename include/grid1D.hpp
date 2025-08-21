#ifndef GRID_1D
#define GRID_1D

#include "function.hpp"
#include "grid.hpp"
#include "cell.hpp"

#include <iostream>
#include <vector>

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

    const double getStart() const {
        return start;
    }
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

    // Returns a 1D grid with boundary tags 0, 1
    operator Grid<1>() const {
        std::vector<Cell<1>> cells(getN()-1);
        std::vector<Point<1>> uniqueNodes(getN());

        for (unsigned int i = 0; i < N; ++i) {
            double current_x = (*this)(i);
            double next_x = (*this)(i + 1);
            uniqueNodes.emplace_back(current_x);
            if (i < N - 1) {
                Cell<1>::NodeVector nodeVector = {Point<1>(current_x), Point<1>(next_x)};
                Cell<1>::NodeIndexes nodeIndexes = {i, i + 1};
                cells.emplace_back(Cell<1>(nodeVector, nodeIndexes));
            }
        }

        std::vector<BoundaryCell<0>> boundaryCells(2);
        boundaryCells.emplace_back(BoundaryCell<0>({getStart()}, {0}, 0));
        boundaryCells.emplace_back(BoundaryCell<0>({getEnd()}, {1}, 1));

        Grid<1> grid(cells, uniqueNodes, boundaryCells);

        return grid;
    }
};

#endif