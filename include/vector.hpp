#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>

class Vector{

    std::vector<double> vector;

    public:

    Vector(int N) : vector(N, 0.0) {};

    inline double& operator[](int index){
        return vector[index];
    };

    Vector operator +(const Vector& other) const;
    Vector operator *(const Vector& other) const;
};

#endif // VECTOR_HPP