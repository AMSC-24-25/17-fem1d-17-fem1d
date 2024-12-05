#include "../include/thomas.hpp"
#include <math.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

typedef SparseMatrix<double, RowMajor> SparseMat;

void Thomas::ForwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs) {
    int n = b.size();

    for(int i=1; i<n; i++){
        double m = a[i] / b[i-1];
        b[i] -= m * c[i-1]; 
        rhs[i] -= m * rhs[i-1];
    }
}
void Thomas::BackwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs) {
    int n = b.size();

    x[n-1] = rhs[n-1]/b[n-1];
    for(int i=n-2; i>=0; i--){
        x[i] = (rhs[i] - (c[i]*x[i+1])) / b[i];
    }
}

VectorXd Thomas::ThomasAlgorithm(SparseMat A, VectorXd& rhs){

    VectorXd x(A.rows());
    VectorXd a(A.rows()-1);
    VectorXd b(A.rows());
    VectorXd c(A.rows()-1);

    MatrixXd mat(A);

    b = mat.diagonal(0); 
    a = mat.diagonal(-1); 
    c = mat.diagonal(1); 

    ForwardSubstitution(a,b,c,x,rhs);
    BackwardSubstitution(a,b,c,x,rhs);
    return x;
}