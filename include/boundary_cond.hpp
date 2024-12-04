#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "function.hpp"

class BoundaryConds {

    private:
    //0 == D, 1 == N
    bool bc1;
    bool bc2;
    Function value1;
    Function value2;

    public:
    BoundaryConds(bool bc1, bool bc2, Function value1, Function value2):
        bc1(bc1), bc2(bc2), value1(value1), value2(value2)
     {};

    bool getBc1() const{ return bc1;}
    bool getBc2() const{ return bc2;}
    Function getValue1() const{ return value1;}
    Function getValue2() const{ return value2;}

};

#endif BOUNDARY_HPP