#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

class Matrix {

    private:
    int N;
    std::vector<std::vector<double>> matrix;

    public:
    Matrix(int N): N(N), matrix(N, std::vector<double>(N, 0.0)) {};
    
    // Add contribution from only 1 segment of mesh to the actual matrix
    void add_contribution(int pos, Matrix small_matrix);

    inline double& operator ()(int i, int j) {
        return matrix[i][j];
    }
    inline int size() const { return N; }

    Matrix& operator +=(Matrix& b);
    Matrix& operator -=(Matrix& b);
    Matrix operator +(Matrix& b) const;
    Matrix operator -(Matrix& b) const;

};


#endif // MATRIX_HPP