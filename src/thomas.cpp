#include <math.h>
#include "../include/matrix.hpp"
#include "../include/vector.hpp"

void ForwardSubstitution(Vector& a, Vector& b, Vector& c, Vector& x, Vector& rhs) {
    int n = b.size();

    for(int i=1; i<n; i++){
        double m = a[i] / b[i-1];
        b[i] -= m * c[i-1]; 
        rhs[i] -= m * c[i-1];
    }
}
void BackwardSubstitution(Vector& a, Vector& b, Vector& c, Vector& x, Vector& rhs) {
    int n = b.size();

    x[n-1] = rhs[n-1]/b[n-1];
    for(int i=n-2; i>=0; i--){
        x[i] = (rhs[i] - (c[i]*x[i+1])) / b[i];
    }
}

Vector ThomasAlgorithm(Matrix A, Vector& rhs){
    Vector x(A.getSize());
    Vector a = A.getDiagonal(-1);
    Vector b = A.getDiagonal(0);
    Vector c = A.getDiagonal(1);
    ForwardSubstitution(a,b,c,x,rhs);
    BackwardSubstitution(a,b,c,x,rhs);
    return x;
}