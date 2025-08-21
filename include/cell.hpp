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

    // Aggiungi metodi mancanti
    const NodeIndexes& getNodeIndexes() const {
        return this->nodeIndices;
    }

    int getBoundaryId() const {
        return boundary_id;
    }
};

// Returns barycentric coordinates of p with respect to triangle cell (2D)
inline std::array<double, 3> barycentricCoordinates(const Cell<2>& cell, const Point<2>& p) {
    const Point<2>& A = cell[0];
    const Point<2>& B = cell[1];
    const Point<2>& C = cell[2];
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

// ============= GEOMETRIA PER BOUNDARY CELLS (EDGE) =============

inline Point<2> mapToGlobalEdge(const BoundaryCell<1>& edge, double xi) {
    // Mappa da coordinata locale xi ∈ [-1, 1] a coordinate globali
    // xi = -1 -> nodo 0, xi = +1 -> nodo 1
    const Point<2>& p1 = edge[0];
    const Point<2>& p2 = edge[1];
    
    double t = (xi + 1.0) / 2.0;  // trasforma da [-1,1] a [0,1]
    double x = (1.0 - t) * p1[0] + t * p2[0];
    double y = (1.0 - t) * p1[1] + t * p2[1];
    
    return Point<2>(x, y);
}

inline Point<2> edgeNormal(const BoundaryCell<1>& edge) {
    // Calcola il vettore normale unitario al segmento (outward normal)
    const Point<2>& p1 = edge[0];
    const Point<2>& p2 = edge[1];
    
    // Vettore tangente
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    
    // Vettore normale (ruotato di 90° in senso orario: (dx,dy) -> (dy,-dx))
    // Questo dà la normale esterna per mesh con orientazione antioraria
    double nx = dy;
    double ny = -dx;
    
    // Normalizza
    double length = sqrt(nx*nx + ny*ny);
    if (length < 1e-14) {
        std::cerr << "Edge degenere in edgeNormal\n";
        exit(-1);
    }
    
    return Point<2>(nx/length, ny/length);
}

inline void getShapeFunctions1D(double xi, std::vector<double>& phi) {
    // Shape functions lineari 1D in [-1, 1]
    phi.resize(2);
    phi[0] = (1.0 - xi) / 2.0;  // shape function per nodo 0
    phi[1] = (1.0 + xi) / 2.0;  // shape function per nodo 1
}


#endif