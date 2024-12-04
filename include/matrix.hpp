#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>

class Matrix {

    private:
    std::vector<std::vector<double>> matrix;
    int size;

    public:
    Matrix(int N): matrix(N, std::vector<double>(N, 0.0)), size(N) {};
    
    // Add contribution from only 1 segment of mesh to the actual matrix
    void add_contributions(Matrix small_matrix, int a);

    inline double& operator ()(int i, int j) {
        return matrix[i][j];
    }

    int getSize() const{ return size;}

};


#endif // MATRIX_HPP