#include "cell.hpp"

template<>
double Cell<1>::measure() const {
    // if (getN() != 2) {
    //     std::cerr << "Cell<1> must be a line\n";
    //     exit(-1);
    // }
    const Cell<1>& cell = *this;
    const Point<1>& A = cell[0];
    const Point<1>& B = cell[1];
    return std::abs(B[0] - A[0]);
}
template<>
double Cell<2>::measure() const {
    if (getN() != 3) {
        std::cerr << "Cell must be a triangle\n";
        exit(-1);
    }

    // Area formula: A = 0.5 * |AB x AC|
    const Cell<2>& cell = *this;
    const Point<2> AB = cell[1] - cell[0];
    const Point<2> AC = cell[2] - cell[0];

    double detT = (AB[0] * AC[1]) - (AC[0] * AB[1]);
    return 0.5 * std::abs(detT);
}
template<>
double Cell<3>::measure() const {
    if (getN() != 4) {
        std::cerr << "Cell must be a tetrahedron\n";
        exit(-1);
    }

    // Volume formula: V = |(AB x AC) · AD| / 6
    const Cell<3>& cell = *this;
    Point<3> AB = cell[1] - cell[0]; 
    Point<3> AC = cell[2] - cell[0]; 
    Point<3> AD = cell[3] - cell[0]; 

    // Cross product AB x AC
    Point<3> cross{
        AB[1] * AC[2] - AB[2] * AC[1],
        AB[2] * AC[0] - AB[0] * AC[2],
        AB[0] * AC[1] - AB[1] * AC[0]
    };

    return std::abs(cross * AD) / 6.0;
}

template<>
Point<1> Cell<1>::barycentricGradient(int i) const {
    if (getN() != 2) {
        std::cerr << "Cell<1> must have exactly 2 nodes\n";
        exit(-1);
    }

    if (i == 0)
        return -1.0 / measure();
    else if (i == 1)
        return 1.0 / measure();

    std::cerr << "Invalid shape function index in barycentricGradient\n";
    exit(-1);
}

// Returns gradient of i-th barycentric shape function (constant on triangles)
template<>
Point<2> Cell<2>::barycentricGradient(int i) const {
    if (getN() != 3) {
        std::cerr << "Cell must be a triangle\n";
        exit(-1);
    }
    const Cell<2>& cell = *this;
    const Point<2>& A = cell[0];
    const Point<2>& B = cell[1];
    const Point<2>& C = cell[2];
    double area2 = measure() * 2.0; // 2 * area

    if (i == 0)
        return Point<2>((B[1] - C[1]) / area2, (C[0] - B[0]) / area2);
    else if (i == 1)
        return Point<2>((C[1] - A[1]) / area2, (A[0] - C[0]) / area2);
    else if (i == 2)
        return Point<2>((A[1] - B[1]) / area2, (B[0] - A[0]) / area2);

    std::cerr << "Indice shape function non valido in barycentricGradient\n";
    exit(-1);
}

// Returns gradient of i-th barycentric shape function (constant on tetrahedra)
template<>
Point<3> Cell<3>::barycentricGradient(int i) const {
    // Gradient of barycentric coordinates in a tetrahedron
    // lambda_i(x) = a_i*x + b_i*y + c_i*z + d_i
    // grad(lambda_i) = (a_i, b_i, c_i)
    if (getN() != 4) {
        std::cerr << "Cell must be a tetrahedron in barycentricGradient\n";
        exit(-1);
    }

    // Tetrahedron determinant
    double detT = measure() * 6.0;

    if (std::abs(detT) < 1e-15) {
        std::cerr << "WARNING: Degenerate tetrahedron (det = "<< std::abs(detT) <<") in barycentricGradient\n";
        return Point<3>(0,0,0);
    }

    // Standard formula: grad(λᵢ) = (face_normal_opposite_to_i) / (6*volume)
    // For tetrahedron ABCD, the gradient of λᵢ is the outward normal to face opposite node i
    const Cell<3>& cell = *this;
    
    // Define the three vertices of the face opposite to node i
    int face_nodes[4][3] = {
        {1, 2, 3},  // Face opposite to node 0: nodes 1,2,3 (BCD)
        {0, 3, 2},  // Face opposite to node 1: nodes 0,3,2 (ADC) - reversed for orientation
        {0, 1, 3},  // Face opposite to node 2: nodes 0,1,3 (ABD)
        {0, 2, 1}   // Face opposite to node 3: nodes 0,2,1 (ACB) - reversed for orientation
    };
    
    const Point<3>& P = cell[face_nodes[i][0]];
    const Point<3>& Q = cell[face_nodes[i][1]];
    const Point<3>& R = cell[face_nodes[i][2]];
    
    // Compute face normal as (Q-P) × (R-P)
    Point<3> PQ = Q - P;
    Point<3> PR = R - P;
    Point<3> normal = Point<3>{
        PQ[1]*PR[2] - PQ[2]*PR[1],
        PQ[2]*PR[0] - PQ[0]*PR[2],
        PQ[0]*PR[1] - PQ[1]*PR[0]
    };
    
    return Point<3>(std::vector<double>{normal[0]/detT, normal[1]/detT, normal[2]/detT});
    std::cerr << "Shape function index " << i << " invalid in barycentricGradient (3D)\n";
    exit(-1);
}

template<>
double BoundaryCell<2>::measure() const {
    // Computes the area of the triangular face
    const Point<3>& p1 = this->getNode(0);
    const Point<3>& p2 = this->getNode(1);
    const Point<3>& p3 = this->getNode(2);

    // Two vectors on the sides of the triangle
    Point<3> v1(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
    Point<3> v2(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
    
    // Norm of the cross product / 2
    double nx = v1[1] * v2[2] - v1[2] * v2[1];
    double ny = v1[2] * v2[0] - v1[0] * v2[2];
    double nz = v1[0] * v2[1] - v1[1] * v2[0];
    
    return 0.5 * sqrt(nx*nx + ny*ny + nz*nz);
}

template<>
double BoundaryCell<1>::measure() const {
    // Computes the length of the edge
    const Point<2>& p1 = this->getNode(0);
    const Point<2>& p2 = this->getNode(1);
    return p2.distance(p1);
}