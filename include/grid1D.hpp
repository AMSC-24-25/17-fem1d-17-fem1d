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
    inline unsigned int getN() const noexcept { return N; }

    // Returns a 1D grid with boundary tags 0, 1
    operator Grid<1>() const {
        // Inizializza vettori vuoti invece di preallocare
        std::vector<Cell<1>> cells;
        std::vector<Point<1>> uniqueNodes;
        
        // Riserva spazio ma non crea elementi
        cells.reserve(getN()-1);
        uniqueNodes.reserve(getN());

        // Aggiungi prima tutti i nodi
        for (unsigned int i = 0; i < N; ++i) {
            double x = (*this)(i);
            uniqueNodes.emplace_back(x);
        }

        // Poi crea celle tra nodi adiacenti
        for (unsigned int i = 0; i < N-1; ++i) {
            Cell<1>::NodeVector nodeVector = {Point<1>((*this)(i)), Point<1>((*this)(i+1))};
            Cell<1>::NodeIndexes nodeIndexes = {i, i+1};
            cells.emplace_back(nodeVector, nodeIndexes);
        }

        // Boundary cells
        std::vector<BoundaryCell<0>> boundaryCells;
        boundaryCells.reserve(2);
        boundaryCells.emplace_back(BoundaryCell<0>({getStart()}, {0}, 0));
        boundaryCells.emplace_back(BoundaryCell<0>({getEnd()}, {getN()-1}, 1));

        return Grid<1>(cells, uniqueNodes, boundaryCells);
    }
};

#endif