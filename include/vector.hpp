#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <stdexcept>
#include <utility>

// Forward declaration to avoid circular include
template <unsigned int dim>
class Function;

template <unsigned int dim>
class Vector {

    private:
    std::vector<Function<dim>> vector;

    public:
    Vector(std::vector<Function<dim>> vec) : vector(std::move(vec)) {};

    size_t size() const { return vector.size(); }

    Vector<dim> operator +(const Vector<dim>& v) const {
        std::vector<Function<dim>> result;
        if (vector.size() != v.vector.size()) {
            throw std::invalid_argument("Vectors must be of the same size");
        }

        for (size_t i = 0; i < vector.size(); ++i) {
            result.push_back(vector[i] + v.vector[i]);
        }
        return Vector<dim>(result);
    }

    Vector<dim> operator -(const Vector<dim>& v) const {
        std::vector<Function<dim>> result;
        if (vector.size() != v.vector.size()) {
            throw std::invalid_argument("Vectors must be of the same size");
        }
        for (size_t i = 0; i < vector.size(); ++i) {
            // If Function::operator- is not defined, emulate as a + (b * -1)
            result.push_back(vector[i] + (v.vector[i] * (-1.0)));
        }
        return Vector<dim>(result);
    }

    Function<dim> operator *(const Vector<dim>& v) const {
        Function<dim> result([](Point<dim> p) { return 0.0; });
        for (size_t i = 0; i < vector.size(); ++i) {
            result += vector[i] * v.vector[i];
        }
        return result;
    }

    Function<dim> operator [](int index) const {
        if (index >= vector.size()) {
            throw std::out_of_range("Index out of range");
        }
        return vector[index];
    }

};

#endif