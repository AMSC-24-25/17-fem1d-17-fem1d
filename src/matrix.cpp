#include "matrix.hpp"

#include <iostream>

void Matrix::add_contribution(int pos, Matrix small_matrix) {

    // TODO: forse da implementare

}

Matrix& Matrix::operator+=(Matrix &b)
{
    if(N != b.size())
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
    if(N != b.size())
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