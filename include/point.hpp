#ifndef POINT
#define POINT

#include <iostream>
#include <vector>

using CoordinateVector = std::vector<double>;

template<unsigned int dim>
struct Point{
    CoordinateVector coords;

    Point(CoordinateVector coords_) : coords(coords_) 
    {}
    Point() : coords(dim, 0.0) {}
    
    Point(double x, double y) : coords{ x, y } {
        static_assert(dim == 2, "Passing 2 coordinates is only valid for 2D points");
    }
    Point(double x, double y, double z) : coords{ x, y, z } {
        static_assert(dim == 3, "Passing 3 coordinates is only valid for 3D points");
    }


    double distance(Point p){
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
};

#endif