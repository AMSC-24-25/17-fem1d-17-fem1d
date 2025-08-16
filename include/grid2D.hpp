#ifndef GRID_2D
#define GRID_2D

#include "cell.hpp"
#include "point.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

class Grid2D{
    using CellVector = std::vector<Cell<2>>;
    using NodeVector = Cell<2>::NodeVector;
    CellVector cells;

public:
    Grid2D(CellVector cells) : cells(cells)
    {}

    Grid2D() : cells() {}

    void addCell(const Cell<2>& cell) {
        cells.push_back(cell);
    }
    unsigned int getNumElements() const {
        return cells.size();
    }
    unsigned int getNumNodes() const {
        return cells.empty() ? 0 : cells[0].getN() * cells.size();
    }
    const Cell<2>& getCell(int i) const {
        if (i < 0 || i >= cells.size()) {
            std::cerr << "Index out of bounds in Grid2D::getCell\n";
            exit(-1);
        }
        return cells[i];
    }

    void parseFromMsh(const std::string& filename);
};

#endif