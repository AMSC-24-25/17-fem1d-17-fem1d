#include "../include/matrix.hpp"

void Matrix::add_contributions(Matrix small_matrix, int a) {
    (*this)(a, a) += small_matrix(0, 0);
    (*this)(a, a+1) += small_matrix(0, 1);
    (*this)(a+1, a) += small_matrix(1, 0);
    (*this)(a+1, a+1) += small_matrix(1, 1);
};