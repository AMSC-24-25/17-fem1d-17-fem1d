#include "cell.hpp"

template<>
Point<1> Cell<1>::barycentricGradient(int i) const {
    if (getN() != 2) {
        std::cerr << "Cell<1> must have exactly 2 nodes\n";
        exit(-1);
    }

    if (i == 0)
        return -1.0;
    else if (i == 1) 
        return 1.0;

    std::cerr << "Indice shape function non valido in barycentricGradient\n";
    exit(-1);
}

// Returns gradient of i-th barycentric shape function (constant on triangles)
template<>
Point<2> Cell<2>::barycentricGradient(int i) const {
    const Cell<2>& cell = *this;
    const Point<2>& AB = cell[1] - cell[0];
    const Point<2>& AC = cell[2] - cell[0];

    double detT = AB[0]*AC[1] - AB[1]*AC[0];
    if (std::abs(detT) < 1e-14) {
        std::cerr << "Triangolo degenere in barycentricGradient\n";
        exit(-1);
    }
    if (i >= 0 && i < 3) {
        // Indici dei nodi opposti a i
        unsigned int idx[2];
        int count = 0;
        for (int j = 0; j < 3; ++j) {
            if (j != i) idx[count++] = j;
        }
        Point<2> oppositeSide = cell[idx[1]] - cell[idx[0]];

        return Point<2>(
            std::vector<double>{
                (oppositeSide[1]) / detT, 
                -(oppositeSide[0]) / detT
            });
    }
    std::cerr << "Indice shape function non valido in barycentricGradient\n";
    exit(-1);
}

// Returns gradient of i-th barycentric shape function (constant on tetrahedra)
template<>
Point<3> Cell<3>::barycentricGradient(int i) const {
    // Gradiente delle coordinate baricentriche in un tetraedro
    // lambda_i(x) = a_i*x + b_i*y + c_i*z + d_i
    // grad(lambda_i) = (a_i, b_i, c_i)
    if (getN() != 4) {
        std::cerr << "Cell must be a tetrahedron in barycentricGradient\n";
        exit(-1);
    }
    const Cell<3>& cell = *this;
    const Point<3>& A = cell[0];
    const Point<3>& B = cell[1];
    const Point<3>& C = cell[2];
    const Point<3>& D = cell[3];

    Point<3> CB = B - C;
    Point<3> DC = C - D;
    Point<3> BD = D - B;

    // Determinante del tetraedro
    double detT = 
        A[0]*(D[1]*CB[2] + B[1]*DC[2] + C[1]*BD[2])
      - A[1]*(B[0]*DC[2] + C[0]*BD[2] + D[0]*CB[2])
      + A[2]*(B[0]*DC[1] + C[0]*BD[1] + D[0]*CB[1]);

    if (std::abs(detT) < 1e-14) {
        std::cerr << "Tetraedro degenere in barycentricGradient\n";
        exit(-1);
    }

    // Gradiente di lambda_i (opposto al nodo i)
    if (i >= 0 && i < 4) {
        // Indici dei nodi opposti a i
        int idx[3];
        int count = 0;
        for (int j = 0; j < 4; ++j) {
            if (j != i) idx[count++] = j;
        }
        const Point<3>& P = cell[idx[0]];
        const Point<3>& Q = cell[idx[1]];
        const Point<3>& R = cell[idx[2]];

        double a = (P[1]*(Q[2]-R[2]) + Q[1]*(R[2]-P[2]) + R[1]*(P[2]-Q[2]));
        double b = (P[0]*(R[2]-Q[2]) + Q[0]*(P[2]-R[2]) + R[0]*(Q[2]-P[2]));
        double c = (P[0]*(Q[1]-R[1]) + Q[0]*(R[1]-P[1]) + R[0]*(P[1]-Q[1]));
        return Point<3>(std::vector<double>{a/detT, b/detT, c/detT});
    }
    std::cerr << "Shape function index " << i << " invalid in barycentricGradient (3D)\n";
    exit(-1);
}

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

    // Formula is A = 0.5 * |AB x AC|
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

    // Formula is V = |(AB x AC) · AD| / 6
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
