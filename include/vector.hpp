#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <stdexcept>
#include <utility>

// forward declaration con due parametri
template <unsigned int dim, unsigned int returnDim>
class Function;

template <unsigned int dim>
class FunctionVector {

private:
    // gradiente: vettore di funzioni scalari (returnDim=1)
    std::vector<Function<dim, 1>> data;

public:
    explicit FunctionVector(std::vector<Function<dim, 1>> vec) : data(std::move(vec)) {}

    size_t size() const { return data.size(); }

    const Function<dim, 1>& operator[](size_t i) const { 
        return data[i]; 
    }
    Function<dim, 1>& operator[](size_t i) { 
        return data[i]; 
    }

    FunctionVector<dim> operator+(const FunctionVector<dim>& v) const {
        if (data.size() != v.data.size()) {
            throw std::invalid_argument("Vectors must be of the same size");
        }
        std::vector<Function<dim, 1>> result;
        result.reserve(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            result.push_back(data[i] + v.data[i]);
        }
        return FunctionVector<dim>(std::move(result));
    }

    FunctionVector<dim> operator-(const FunctionVector<dim>& v) const {
        if (data.size() != v.data.size()) {
            throw std::invalid_argument("Vectors must be of the same size");
        }
        std::vector<Function<dim, 1>> result;
        result.reserve(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            result.push_back(data[i] + (v.data[i] * -1.0));
        }
        return FunctionVector<dim>(std::move(result));
    }
};

#endif