#ifndef CELL
#define CELL

#include <iostream>
#include <vector>

template<unsigned int dim>
struct Cell {
    using NodeVector = std::vector<Point<dim>>;
    NodeVector nodes;

    Cell(NodeVector nodes_) : nodes(nodes_)
    {}

    unsigned int getN() const {
        return nodes.size();
    }
    const Point<dim>& operator[](unsigned int i) const {
        if (i >= getN()) {
            std::cerr << "Index out of bounds in Cell\n";
            exit(-1);
        }
        return nodes[i];
    }
};

#endif