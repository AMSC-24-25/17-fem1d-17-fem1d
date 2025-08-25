/**
 * @file point.hpp
 * @brief N-dimensional point class for geometric computations
 */
#ifndef POINT
#define POINT

#include <iostream>
#include <cmath>
#include <vector>
#include <array>

using CoordinateVector = std::vector<double>;

/**
 * @brief N-dimensional point with coordinate access and arithmetic operations
 */
template<unsigned int dim>
class Point{

protected: CoordinateVector coords;

public:

    static Point<dim> zero() {
        return Point<dim>(std::vector<double>(dim,0.0));
    }
    static Point<dim> one() {
        return Point<dim>(std::vector<double>(dim,1.0));
    }


    Point(CoordinateVector coords_) : coords(coords_) 
    {}
    Point(std::array<double, dim> coords_) : coords(coords_.begin(), coords_.end()) 
    {}
    Point() : coords(dim, 0.0) {}
    
    Point(double x) : coords{ x } {
        static_assert(dim == 1, "Passing 1 coordinate is only valid for 1D points");
    }
    Point(double x, double y) : coords{ x, y } {
        static_assert(dim == 2, "Passing 2 coordinates is only valid for 2D points");
    }
    Point(double x, double y, double z) : coords{ x, y, z } {
        static_assert(dim == 3, "Passing 3 coordinates is only valid for 3D points");
    }

    operator double() const {
        static_assert(dim == 1, "Automatic cast to double is only valid for Point<1>");
        return coords[0];
    }

    double distance(Point p) const {
        double sum = 0.0;
        for (unsigned int i = 0; i < dim; ++i) {
            double diff = coords[i] - p.coords[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }

    double operator[](unsigned int i) const {
        if (i >= dim) {
            std::cerr << "Index out of bounds in Point\n";
            exit(-1);
        }
        return coords[i];
    }
    Point<dim> operator-(Point<dim> other) const {
        CoordinateVector result(dim);
        for (unsigned int i = 0; i < dim; ++i) {
            result[i] = coords[i] - other[i];
        }
        return Point<dim>(result);
    }
    Point<dim> operator+(Point<dim> other) const {
        CoordinateVector result(dim);
        for (unsigned int i = 0; i < dim; ++i) {
            result[i] = coords[i] + other[i];
        }
        return Point<dim>(result);
    }

    // Dot product
    double operator *(Point<dim> p) const {
        double result = 0.0;
        for (unsigned int i = 0; i < dim; ++i) {
            result += coords[i] * p.coords[i];
        }
        return result;
    }
    Point<dim> operator *(double scalar) const {
        Point<dim> result(coords);
        for (unsigned int i = 0; i < dim; ++i) {
            result.coords[i] *= scalar;
        }
        return result;
    }

    double x() const {
        static_assert(dim >= 1, "Point is not 1D");
        return coords[0];
    }

    double y() const {
        static_assert(dim >= 2, "Point is not 2D");
        return coords[1];
    }

    double z() const {
        static_assert(dim >= 3, "Point is not 3D");
        return coords[2];
    }

};

// Product between points of different dimensions
template<unsigned int dim1, unsigned int dim2>
Point< (dim1 == 1 ? dim2 : (dim2 == 1 ? dim1 : 1)) > operator*(const Point<dim1>& p1, const Point<dim2>& p2) {
    if constexpr (dim1 == 1) {
    // Scalar-vector product: returns a vector of the same dimension as the second operand
        std::array<double, dim2> coords;
        for (int i = 0; i < dim2; i++)
            coords[i] = p1.x() * p2[i];
        return Point<dim2>(coords);
    } 
    else if constexpr (dim2 == 1) {
    // Vector-scalar product: returns a vector of the same dimension as the first operand
        std::array<double, dim1> coords;
        for (int i = 0; i < dim1; i++)
            coords[i] = p1[i] * p2.x();
        return Point<dim1>(coords);
    } 
    else if constexpr (dim1 == dim2) {
    // Dot product between vectors of the same dimension: returns a scalar
        double result = 0.0;
        for (unsigned int i = 0; i < dim1; ++i) {
            result += p1[i] * p2[i];
        }
        return Point<1>(result);
    }
    else {
    // For completeness, this branch should never be reached
        static_assert(dim1 == dim2, "Incompatible dimensions for multiplication");
    return Point<1>(); // To avoid warning "not all code paths return a value"
    }
}

#endif