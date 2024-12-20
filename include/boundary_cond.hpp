#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "function.hpp"


class BoundaryCond {

    private:
    bool is_neumann; //False == D, True == N
    Function boundary;

    public:
    BoundaryCond(bool is_neumann, Function boundary) : is_neumann(is_neumann), boundary(boundary) {};

    inline const bool isNeumann() const{
        return is_neumann;
    }

    inline const Function getBoundary() const{
        return boundary;
    }

};

#endif