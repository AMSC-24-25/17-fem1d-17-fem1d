#include "../include/matrix.hpp"

#include <iostream>

void Matrix::add_contribution(int a, Matrix small_matrix) {
    (*this)(a, a) += small_matrix(0, 0);
    (*this)(a, a+1) += small_matrix(0, 1);
    (*this)(a+1, a) += small_matrix(1, 0);
    (*this)(a+1, a+1) += small_matrix(1, 1);
}

Matrix& Matrix::operator+=(Matrix &b)
{
    if(N != b.getSize())
        std::cerr << "b.size != this.size";


    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            matrix[i][j] += b.matrix[i][j];
        }

    return *this;
}

Matrix& Matrix::operator-=(Matrix &b)
{
    if(N != b.getSize())
        std::cerr << "b.size != this.size";

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            if(b(i,j) != 0)
                (*this)(i,j) += b(i,j);
        }

    return *this;
}

Matrix Matrix::operator +(Matrix& b) const {
    Matrix result(*this);
    result += b;
    return result;
}
Matrix Matrix::operator -(Matrix& b) const {
    Matrix result(*this);
    result -= b;
    return result;
}

Vector Matrix::getDiagonal(int diag) const {
    if (N == 0) {
        std::cerr << "Empty matrix";
    }

    int diagonalSize = 0;
    if (diag >= 0) {
        for (int i = 0; i < N && i + diag < N; ++i) {
            ++diagonalSize;
        }
    } else {
        for (int i = 0; i < N && i - diag < N; ++i) {
            ++diagonalSize;
        }
    }

    Vector diagonal(diagonalSize);
    for (int i = 0; i < diagonalSize; ++i) {
        int row, col;
        if (diag >= 0) {
            row = i;
            col = i + diag;
        } else {
            row = i - diag;
            col = i;
        }
        diagonal[i] = (*this)(row,col);
    }

    return diagonal;
}