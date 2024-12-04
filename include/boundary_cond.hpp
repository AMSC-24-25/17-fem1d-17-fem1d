#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "function.hpp"

class BoundaryConds {

    private: 
    bool bc1;
    bool bc2;
    Function value1;
    Function value2;

    public:
    BoundaryConds(bool bc1, bool bc2, Function value1, Function value2):
        bc1(bc1), bc2(bc2), value1(value1), value2(value2)
     {};

};

#endif