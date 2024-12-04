#ifndef THOMAS
#define THOMAS

#include "vector.hpp"
#include "matrix.hpp"

class Thomas {
    
    public:

    Thomas() {};
    
    void ForwardSubstitution(Vector& a, Vector& b, Vector& c, Vector& x, Vector& rhs);

    void BackwardSubstitution(Vector& a, Vector& b, Vector& c, Vector& x, Vector& rhs);
    
    Vector ThomasAlgorithm(Matrix A, Vector& rhs);

};

#endif // THOMAS