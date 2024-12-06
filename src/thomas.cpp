#include "../include/thomas.hpp"

using namespace Eigen;

typedef SparseMatrix<double, RowMajor> SparseMat;

void Thomas::ForwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs) {
    int n = b.size();

    for(int i=1; i<n; i++){
        if(b[i-1] != 0) throw std::runtime_error("A diag has zero value");
        
        double m = a[i-1] / b[i-1];
        b[i] -= m * c[i-1]; 
        rhs[i] -= m * rhs[i-1];
    }

}


void Thomas::BackwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs) {
    int n = b.size();

    // std::cout << "B = " << b[n-1] << std::endl;
    // std::cout << "RHS = " << b[n-1] << std::endl;
    x[n-1] = b[n-1] != 0 ? rhs[n-1]/b[n-1] : throw std::runtime_error("rhs has zero value");
    // std::cout << "X = " << x[n-1] << std::endl;
    for(int i=n-2; i>=0; i--){
        // std::cout << "B[" << i << "] = " << i << b[i] << std::endl;
        x[i] = b[i] != 0 ? (rhs[i] - (c[i]*x[i+1])) / b[i] : throw std::runtime_error("rhs has zero value");
    }
}

VectorXd Thomas::ThomasAlgorithm(SparseMat A, VectorXd& rhs){

    VectorXd x = VectorXd::Zero(A.rows());
    VectorXd a(A.rows()-1);
    VectorXd b(A.rows());
    VectorXd c(A.rows()-1);

    std::unique_ptr<MatrixXd> mat = std::make_unique<MatrixXd>(MatrixXd(A));

    b = VectorXd(mat->diagonal(0)); 
    a = VectorXd(mat->diagonal(-1)); 
    c = VectorXd(mat->diagonal(1)); 

    mat.reset();

    ForwardSubstitution(a,b,c,x,rhs);
    BackwardSubstitution(a,b,c,x,rhs);
    // std::cout << "X:\n" << x << std::endl;
    return x;
}