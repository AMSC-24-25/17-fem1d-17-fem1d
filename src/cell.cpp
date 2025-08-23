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
    const Point<3>& A = cell[0];
    const Point<3>& B = cell[1]; 
    const Point<3>& C = cell[2];
    const Point<3>& D = cell[3];

    if (i == 0) {
        // grad(λ₀): normal to face BCD (opposite to A)
        Point<3> BC = C - B;
        Point<3> BD = D - B;
        Point<3> normal = Point<3>{
            BC[1]*BD[2] - BC[2]*BD[1],
            BC[2]*BD[0] - BC[0]*BD[2],
            BC[0]*BD[1] - BC[1]*BD[0]
        };
        return Point<3>(std::vector<double>{normal[0]/detT, normal[1]/detT, normal[2]/detT});
    } else if (i == 1) {
        // grad(λ₁): normal to face ACD (opposite to B)
        Point<3> AC = C - A;
        Point<3> AD = D - A;
        Point<3> normal = Point<3>{
            AD[1]*AC[2] - AD[2]*AC[1],  // Note: AD × AC (reversed for correct orientation)
            AD[2]*AC[0] - AD[0]*AC[2],
            AD[0]*AC[1] - AD[1]*AC[0]
        };
        return Point<3>(std::vector<double>{normal[0]/detT, normal[1]/detT, normal[2]/detT});
    } else if (i == 2) {
        // grad(λ₂): normal to face ABD (opposite to C)
        Point<3> AB = B - A;
        Point<3> AD = D - A;
        Point<3> normal = Point<3>{
            AB[1]*AD[2] - AB[2]*AD[1],
            AB[2]*AD[0] - AB[0]*AD[2],
            AB[0]*AD[1] - AB[1]*AD[0]
        };
        return Point<3>(std::vector<double>{normal[0]/detT, normal[1]/detT, normal[2]/detT});
    } else if (i == 3) {
        // grad(λ₃): normal to face ABC (opposite to D)
        Point<3> AB = B - A;
        Point<3> AC = C - A;
        Point<3> normal = Point<3>{
            AC[1]*AB[2] - AC[2]*AB[1],  // Note: AC × AB (reversed for correct orientation)
            AC[2]*AB[0] - AC[0]*AB[2],
            AC[0]*AB[1] - AC[1]*AB[0]
        };
        return Point<3>(std::vector<double>{normal[0]/detT, normal[1]/detT, normal[2]/detT});
    }
    std::cerr << "Shape function index " << i << " invalid in barycentricGradient (3D)\n";
    exit(-1);
}
