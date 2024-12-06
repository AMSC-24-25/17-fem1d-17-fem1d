#ifndef THOMAS
#define THOMAS

#include <math.h>
#include <memory>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

typedef SparseMatrix<double, RowMajor> SparseMat;

class Thomas {
    
    public:

    Thomas() {};
    
    void ForwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs);

    void BackwardSubstitution(VectorXd& a, VectorXd& b, VectorXd& c, VectorXd& x, VectorXd& rhs);
    
    VectorXd ThomasAlgorithm(SparseMat A, VectorXd& rhs);

};

#endif // THOMAS