#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

class Matrix {

    private:
    std::vector<std::vector<double>> matrix;

    public:
    Matrix(int N): matrix(N, std::vector<double>(N, 0.0)) {};
    
    // Add contribution from only 1 segment of mesh to the actual matrix
    void add_contributions(Matrix small_matrix);

};


#endif // MATRIX_HPP