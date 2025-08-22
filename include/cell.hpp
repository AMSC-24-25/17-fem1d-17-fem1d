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
        std::cerr << "Degenerate triangle in barycentricCoordinates\n";
        exit(-1);
    }
    double l1 = ((xB - x)*(yC - y) - (xC - x)*(yB - y)) / detT;
    double l2 = ((xC - x)*(yA - y) - (xA - x)*(yC - y)) / detT;
    double l3 = 1.0 - l1 - l2;
    return {l1, l2, l3};
}

// ============= GEOMETRY FOR BOUNDARY CELLS (EDGE) =============

inline Point<2> mapToGlobalEdge(const BoundaryCell<1>& edge, double xi) {
    // Maps from local coordinate xi ∈ [0, 1] to global coordinates
    const Point<2>& p1 = edge[0];
    const Point<2>& p2 = edge[1];
    double x = (1.0 - xi) * p1[0] + xi * p2[0];
    double y = (1.0 - xi) * p1[1] + xi * p2[1];
    return Point<2>(x, y);
}

inline Point<2> edgeNormal(const BoundaryCell<1>& edge) {
    // Computes the unit normal vector to the segment (outward normal)
    const Point<2>& p1 = edge[0];
    const Point<2>& p2 = edge[1];
    
    // Tangent vector
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    
    // Normal vector (rotated 90° clockwise: (dx,dy) -> (dy,-dx))
    // This gives the outward normal for meshes with counterclockwise orientation
    double nx = dy;
    double ny = -dx;
    
    // Normalize
    double length = sqrt(nx*nx + ny*ny);
    if (length < 1e-14) {
        std::cerr << "Degenerate edge in edgeNormal\n";
        exit(-1);
    }
    
    return Point<2>(nx/length, ny/length);
}

inline void getShapeFunctions1D(double xi, std::vector<double>& phi) {
    // Linear 1D shape functions in [0, 1]
    phi.resize(2);
    phi[0] = 1.0 - xi;  // shape function for node 0
    phi[1] = xi;        // shape function for node 1
}

// ============= GEOMETRY FOR BOUNDARY CELLS (3D FACES) =============

inline Point<3> mapToGlobalFace(const BoundaryCell<2>& face, double xi, double eta) {
    // Maps from local coordinates (xi, eta) to global coordinates for a triangular face
    // Barycentric coordinate: (1-xi-eta, xi, eta)
    const Point<3>& p1 = face[0];
    const Point<3>& p2 = face[1];
    const Point<3>& p3 = face[2];
    
    double x = (1.0 - xi - eta) * p1[0] + xi * p2[0] + eta * p3[0];
    double y = (1.0 - xi - eta) * p1[1] + xi * p2[1] + eta * p3[1];
    double z = (1.0 - xi - eta) * p1[2] + xi * p2[2] + eta * p3[2];
    
    return Point<3>(x, y, z);
}

inline Point<3> faceNormal(const BoundaryCell<2>& face) {
    // Computes the unit normal vector to the triangular face
    const Point<3>& p1 = face[0];
    const Point<3>& p2 = face[1];
    const Point<3>& p3 = face[2];
    
    // Two vectors on the sides of the triangle
    Point<3> v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
    Point<3> v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
    
    // Cross product v1 × v2
    double nx = v1[1] * v2[2] - v1[2] * v2[1];
    double ny = v1[2] * v2[0] - v1[0] * v2[2];
    double nz = v1[0] * v2[1] - v1[1] * v2[0];
    
    // Normalize
    double length = sqrt(nx*nx + ny*ny + nz*nz);
    if (length < 1e-14) {
        std::cerr << "Degenerate face in faceNormal\n";
        exit(-1);
    }
    
    return Point<3>(nx/length, ny/length, nz/length);
}

inline double faceArea(const BoundaryCell<2>& face) {
    // Computes the area of the triangular face
    const Point<3>& p1 = face[0];
    const Point<3>& p2 = face[1];
    const Point<3>& p3 = face[2];
    
    // Two vectors on the sides of the triangle
    Point<3> v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
    Point<3> v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
    
    // Norm of the cross product / 2
    double nx = v1[1] * v2[2] - v1[2] * v2[1];
    double ny = v1[2] * v2[0] - v1[0] * v2[2];
    double nz = v1[0] * v2[1] - v1[1] * v2[0];
    
    return 0.5 * sqrt(nx*nx + ny*ny + nz*nz);
}

inline void getShapeFunctions2D(double xi, double eta, std::vector<double>& phi) {
    // Linear 2D shape functions for triangle in barycentric coordinates
    phi.resize(3);
    phi[0] = 1.0 - xi - eta;  // shape function for node 0
    phi[1] = xi;              // shape function for node 1
    phi[2] = eta;             // shape function for node 2
}


#endif