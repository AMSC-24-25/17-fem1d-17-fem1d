#include "../include/thomas.hpp"

using namespace Eigen;

typedef SparseMatrix<double, RowMajor> SparseMat;

void Thomas::forwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs) {
    int n = b.size();

    for(int i=1; i<n; i++){
        if(b[i-1] == 0) throw std::runtime_error("A diag has zero value");
        
        double m = a[i-1] / b[i-1];
        b[i] -= m * c[i-1]; 
        rhs[i] -= m * rhs[i-1];
    }

}


void Thomas::backwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs) {
    int n = b.size();
    
    if(b[n-1] == 0) throw std::runtime_error("rhs has zero value");
    x[n-1] = rhs[n-1]/b[n-1];
    
    for(int i=n-2; i>=0; i--){
        if(b[i] == 0) throw std::runtime_error("rhs has zero value");
        x[i] = (rhs[i] - (c[i]*x[i+1])) / b[i];
    }
}


VectorXd Thomas::solve(SparseMat A, VectorXd& rhs){

    VectorXd x = VectorXd::Zero(A.rows());
    VectorXd a(A.rows()-1);
    VectorXd b(A.rows());
    VectorXd c(A.rows()-1);

    std::unique_ptr<MatrixXd> mat = std::make_unique<MatrixXd>(MatrixXd(A));

    b = VectorXd(mat->diagonal(0)); 
    a = VectorXd(mat->diagonal(-1)); 
    c = VectorXd(mat->diagonal(1)); 

    mat.reset();

    forwardSubstitution(a,b,c,x,rhs);
    backwardSubstitution(a,b,c,x,rhs);
    // std::cout << "X:\n" << x << std::endl;
    return x;
}