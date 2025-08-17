#ifndef FUNCTION
#define FUNCTION

#include <functional>
#include <iostream>
#include "point.hpp"


template<unsigned int dim>
class Function{
    using fun = std::function<double(Point<dim>)>;

    static const fun zeroFun = [](Point<dim> p) -> double { return 0; };
    static const fun oneFun =  [](Point<dim> p) -> double { return 1; };

    private:
    const fun function;
    const std::vector<fun> gradient;

    public:
    explicit Function(fun f) : function(f) , gradient{zeroFun} {}

    explicit Function(fun f, fun gx) : function(f), gradient{gx} {
        static_assert(dim >= 1, "Error in gradient: Dimension is less than 1D");
    }
    explicit Function(fun f, fun gx, fun gy) : function(f), gradient{gx, gy} {
        static_assert(dim >= 2, "Error in gradient: Dimension is less than 2D");
    }

    explicit Function(fun f, fun gx, fun gy, fun gz) : function(f), gradient{gx, gy, gz} {
        static_assert(dim >= 3, "Error in gradient: Dimension is less than 3D");
    }


    inline virtual double value(Point<dim> p) const{
        return function(p);
    }
    inline virtual double dx(Point<dim> p) const{
        return gradient[0](p);
    }
    inline virtual double dy(Point<dim> p) const{
        return gradient[1](p);
    }
    inline virtual double dz(Point<dim> p) const{
        return gradient[2](p);
    }

    inline virtual std::vector<fun> gradient(Point<dim> p) const{}


    // vecchio metodo

    // /**
    //  * @return this a new Function evaluating to its gradient
    //  */
    // virtual Function getGrad() const {
    //     return Function(gradient);
    // }

    Function operator +(const Function& f) const;
    Function operator *(const Function& f) const;
    Function operator +(const double k) const;
    Function operator *(const double k) const;
    inline virtual double operator ()(double x) const{
        return value(x);
    }

 };

 class ZeroFunction : public Function {
    public:
    ZeroFunction() : Function(zeroFun, zeroFun) {}
 };
 
 class OneFunction : public Function {
    public:
    OneFunction() : Function(oneFun, zeroFun) {}
 };


#endif //FUNCTION