#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>

class Vector{
    private:

    std::vector<double> vector;

    public:

    Vector(int N) : vector(N, 0.0) {};

    inline const double& operator[](int index) const {
        return vector[index];
    }
    inline double& operator[](int index) {
        return vector[index];
    }

    inline double size() const { return vector.size(); }

    Vector operator +(const Vector& other) const;
    Vector operator *(const Vector& other) const;
};

#endif // VECTOR_HPP