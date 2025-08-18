#ifndef CELL
#define CELL

#include "point.hpp"
#include <iostream>
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
    const Point<dim>& getCell(unsigned int i) const{
        return *this[i];
    }
    const unsigned int& getNodeIndex(unsigned int i) const {
        if (i >= getN()) {
            std::cerr << "Index out of bounds for node index in Cell\n";
            exit(-1);
        }
        return nodeIndices[i];
    }
};

template<unsigned int dim>
struct BoundaryCell : public Cell<dim+1> {
    using NodeVector = typename Cell<dim+1>::NodeVector;
    using NodeIndexes = typename Cell<dim+1>::NodeIndexes;

    unsigned int boundary_id;

    BoundaryCell(NodeVector nodes_, NodeIndexes nodeIndices_, unsigned int boundary_id_) 
    : Cell<dim+1>(nodes_, nodeIndices_), boundary_id(boundary_id_)
    {
        if (nodes_.size() < 2) {
            std::cerr << "BoundaryCell must have at least 2 nodes\n";
            exit(-1);
        }
    }
    BoundaryCell() = default;
};

#include <array>

// Calcola le coordinate baricentriche di p rispetto al triangolo cell (2D)
inline std::array<double, 3> barycentricCoordinates(const Cell<2>& cell, const Point<2>& p) {
    const auto& A = cell[0];
    const auto& B = cell[1];
    const auto& C = cell[2];
    double x = p[0], y = p[1];
    double xA = A[0], yA = A[1];
    double xB = B[0], yB = B[1];
    double xC = C[0], yC = C[1];

    double detT = (xB - xA)*(yC - yA) - (xC - xA)*(yB - yA);
    if (std::abs(detT) < 1e-14) {
        std::cerr << "Triangolo degenere in barycentricCoordinates\n";
        exit(-1);
    }
    double l1 = ((xB - x)*(yC - y) - (xC - x)*(yB - y)) / detT;
    double l2 = ((xC - x)*(yA - y) - (xA - x)*(yC - y)) / detT;
    double l3 = 1.0 - l1 - l2;
    return {l1, l2, l3};
}

// Restituisce il gradiente della shape function baricentrica i-esima (costante su triangolo)
inline Point<2> barycentricGradient(const Cell<2>& cell, int i) {
    const auto& A = cell[0];
    const auto& B = cell[1];
    const auto& C = cell[2];
    double xA = A[0], yA = A[1];
    double xB = B[0], yB = B[1];
    double xC = C[0], yC = C[1];
    double detT = (xB - xA)*(yC - yA) - (xC - xA)*(yB - yA);
    if (std::abs(detT) < 1e-14) {
        std::cerr << "Triangolo degenere in barycentricGradient\n";
        exit(-1);
    }
    // Gradiente di lambda_1
    if (i == 0) {
        return Point<2>(std::vector<double>{(yB - yC)/detT, (xC - xB)/detT});
    }
    // Gradiente di lambda_2
    if (i == 1) {
        return Point<2>(std::vector<double>{(yC - yA)/detT, (xA - xC)/detT});
    }
    // Gradiente di lambda_3
    if (i == 2) {
        return Point<2>(std::vector<double>{(yA - yB)/detT, (xB - xA)/detT});
    }
    std::cerr << "Indice shape function non valido in barycentricGradient\n";
    exit(-1);
}

#endif