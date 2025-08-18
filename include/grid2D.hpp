#ifndef GRID_2D
#define GRID_2D

#include "cell.hpp"
#include "point.hpp"
#include "phi_function2d.hpp"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

class Grid2D{
    using CellVector = std::vector<Cell<2>>;
    using BoundaryCellVector = std::vector<BoundaryCell<1>>;
    using NodeVector = Cell<2>::NodeVector;
    using BoundaryNodeVector = Cell<1>::NodeVector;
    using PhiFunctionVector2D = std::vector<PhiFunction2D>;
    // using BoundaryPhiFunctionVector2D = std::vector<BoundaryPhiFunction2D>;

    CellVector cells;
    BoundaryCellVector boundary_cells;
    NodeVector unique_nodes;

public:
    Grid2D(CellVector cells, NodeVector unique_nodes) : cells(cells), unique_nodes(unique_nodes)
    {}

    Grid2D() : cells(), unique_nodes() {}

    void addCell(const Cell<2>& cell) {
        cells.push_back(cell);
    }
    unsigned int getNumElements() const {
        return cells.size();
    }
    unsigned int getNumNodes() const {
        return unique_nodes.size();
    }
    const Cell<2>& getCell(unsigned int i) const {
        if (i >= cells.size()) {
            std::cerr << "Index out of bounds in Grid2D::getCell\n";
            exit(-1);
        }
        return cells[i];
    }
    const Point<2>& getNode(unsigned int i) const {
        if (i >= unique_nodes.size()) {
            std::cerr << "Index out of bounds in Grid2D::getNode\n";
            exit(-1);
        }
        return unique_nodes[i];
    }

    void parseFromMsh(const std::string& filename);
    
    PhiFunctionVector2D getPhiFunctions() const;
    // PhiFunctionVector2D getBoundaryPhiFunctions() const;
};

#endif