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
        std::cerr << "Not yet implemented Grid2D::getNumNodes." << std::endl;
        // Pay attention: some nodes are shared between cells, so this is not a simple count.
        exit(-1);
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