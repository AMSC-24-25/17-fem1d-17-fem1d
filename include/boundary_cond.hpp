#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "function.hpp"

//Ora va solo per bordi 1D

class BoundaryCond {

    private:
    bool is_neumann; //False == D, True == N
    Function<1> boundary;

    public:
    BoundaryCond(bool is_neumann, Function<1> boundary) : is_neumann(is_neumann), boundary(boundary) {};

    inline const bool isNeumann() const{
        return is_neumann;
    }

    inline const Function<1>& getBoundary() const{
        return boundary;
    }

};

#endif  