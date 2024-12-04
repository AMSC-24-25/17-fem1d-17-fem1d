#include <iostream>
#include "../include/vector.hpp"

Vector Vector::operator +(const Vector& other) const {
    
    if((*this).size() != other.size())
        std::cerr << "Different sizes";

    Vector result(vector.size());
    for (int i = 0; i < size(); ++i) {
        result[i] = (*this)[i] + other[i];
    }

    return result;
}

Vector Vector::operator *(const Vector& other) const{

    if((*this).size() != other.size())
        std::cerr << "Different sizes";

    Vector result(vector.size());
    for (int i = 0; i < size(); ++i) {
        result[i] = (*this)[i] * other[i];
    }

    return result;
}