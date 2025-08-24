#ifndef CELL
#define CELL

#include "point.hpp"
#include <iostream>
#include <array>
#include <vector>

template<unsigned int dim>
struct Cell {
    using NodeVector = std::vector<Point<dim>>;
    using NodeIndexes = std::vector<unsigned int>;
    
    const NodeVector nodes;
    const NodeIndexes nodeIndices;

    Cell(NodeVector nodes_, NodeIndexes nodeIndices_) : nodes(nodes_), nodeIndices(nodeIndices_)
    {
        if (nodes.size() != nodeIndices.size()) {
            std::cerr << "Inconsistent node and index sizes in Cell\n";
            exit(-1);
        }
    }
    Cell() = default;

    unsigned int getN() const {
        return nodes.size();
    }
    const Point<dim>& operator[](unsigned int i) const {
        if (i >= getN()) {
            std::cerr << "Index out of bounds for node in Cell\n";
            exit(-1);
        }
        return nodes[i];
    }
    const Point<dim>& getNode(unsigned int i) const{
        return (*this)[i];
    }
    const unsigned int& getNodeIndex(unsigned int i) const {
        if (i >= getN()) {
            std::cerr << "Index out of bounds for node index in Cell\n";
            exit(-1);
        }
        return nodeIndices[i];
    }
    const NodeIndexes& getNodeIndexes() const {
        return nodeIndices;
    }

    inline Point<dim> mapToGlobal(const Point<dim+1>& barycentricCoords) const {
        const Cell<dim>& cell = *this;
        std::array<double, dim> pCoords;
        for (size_t i = 0; i < dim; i++){
            pCoords[i] = 0;
            for (size_t j = 0; j < dim + 1; j++)
                pCoords[i] += barycentricCoords[j] * cell[j][i];
        }
        return Point<dim>(pCoords);
    }
    
    // Returns gradient of i-th barycentric shape function (constant on tetrahedra)
    Point<dim> barycentricGradient(int i) const;
    double measure() const;
};

template<unsigned int dim>
struct BoundaryCell : public Cell<dim+1> {
    using NodeVector = typename Cell<dim+1>::NodeVector;
    using NodeIndexes = typename Cell<dim+1>::NodeIndexes;

    unsigned int boundary_id;

    BoundaryCell(NodeVector nodes_, NodeIndexes nodeIndices_, unsigned int boundary_id_) 
    : Cell<dim+1>(nodes_, nodeIndices_), boundary_id(boundary_id_) {};
    BoundaryCell() = default;

    unsigned int getBoundaryId() const {
        return boundary_id;
    }

    inline Point<dim+1> mapToGlobal(const Point<dim+1>& barycentricCoords) const {
        const BoundaryCell<dim>& cell = *this;
        std::array<double, dim+1> result;
        for (unsigned int i = 0; i < dim+1; ++i) {
            result[i] = 0.0;
            for (unsigned int j = 0; j < dim+1; ++j) {
                result[i] += barycentricCoords[j] * cell[j][i];
            }
        }
        return Point<dim+1>(result);
    }

    double measure() const;
};

#endif