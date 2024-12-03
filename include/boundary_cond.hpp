#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "function.hpp"


class BoundaryCond {

    private:
    bool is_neuman;
    Function boundary;
    //0 == D, 1 == N

    public:
    BoundaryCond(bool is_neuman, Function boundary) : boundary(boundary), is_neuman(is_neuman) {};

    inline const bool isNeuman() const{
        return is_neuman;
    }

    inline const Function getBoundary() const{
        return boundary;
    }

};

#endif