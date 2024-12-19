#ifndef THOMAS
#define THOMAS

#include <math.h>
#include <memory>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>


using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;

class Thomas {
    
    public:

    Thomas() {};

private:
    
    void forwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs);

    void backwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs);

public:

    VectorXd solve(SparseMat A, VectorXd& rhs);

};

#endif // THOMAS