#ifndef FUNCTION
#define FUNCTION

#include <functional>
#include <iostream>
#include "point.hpp"



template<unsigned int dim> 
class Function{

    using fun = std::function<double(Point<dim>)>;


    static const inline fun zeroFun = [](Point<dim> p) -> double { return 0; };
    static const inline fun oneFun =  [](Point<dim> p) -> double { return 1; };

    private:
    const fun function;
    const std::vector<fun> gradient;

    public:
    explicit Function(fun f) : function(f) , gradient{zeroFun} {}

    explicit Function(fun f, fun gx) : function(f), gradient{gx} {
        static_assert(dim == 1, "Error in gradient: Dimension is less than 1D");
    }
    explicit Function(fun f, fun gx, fun gy) : function(f), gradient{gx, gy} {
        static_assert(dim == 2, "Error in gradient: Dimension is less than 2D");
    }

    explicit Function(fun f, fun gx, fun gy, fun gz) : function(f), gradient{gx, gy, gz} {
        static_assert(dim == 3, "Error in gradient: Dimension is less than 3D");
    }

    explicit Function(fun f, std::vector<fun> gradient) : function(f), gradient(gradient) {
        static_assert(dim == 3, "Error in gradient: Dimension is less than 3D");
    }


    inline double value(Point<dim> p) const{
        return function(p);
    }

    inline double dx_value(Point<dim> p) const{
        return gradient[0](p);
    }
    inline double dy_value(Point<dim> p) const{
        return gradient[1](p);
    }
    inline double dz_value(Point<dim> p) const{
        return gradient[2](p);
    }

    inline std::vector<double> getGradValues(Point<dim> p) const{
        std::vector<double> gradValues;
        for (const fun grad : gradient) {
            gradValues.push_back(grad(p));
        }
        return gradValues;
    }

    inline  std::vector<Function<dim>> getGrad() const{
        std::vector<Function<dim>> gradFunctions;
        for (const fun& grad : gradient) {
            gradFunctions.emplace_back(Function(grad));
        }
        return gradFunctions;
    }


    // vecchio metodo

    // /**
    //  * @return this a new Function evaluating to its gradient
    //  */
    // virtual Function getGrad() const {
    //     return Function(gradient);
    // }

    Function<dim> operator +(const Function<dim>& f) const;
    Function<dim> operator *(const Function<dim>& f) const;
    Function<dim> operator +(const double k) const;
    Function<dim> operator *(const double k) const;

    inline double operator ()(Point<dim> p) const{
        return value(p);
    }

};

#include "function.tpp"


template<unsigned int dim>
class ZeroFunction : public Function<dim> {

    using fun = std::function<double(Point<dim>)>;

    public:
    static const inline fun zeroFun = [](Point<dim> p) -> double { return 0; };


    ZeroFunction() : Function<dim>(zeroFun, zeroFun) {}
 };


template<unsigned int dim>
class OneFunction : public Function<dim> {

    using fun = std::function<double(Point<dim>)>;

    public:
    static const inline fun oneFun = [](Point<dim> p) -> double { return 1; };
    static const inline fun zeroFun = [](Point<dim> p) -> double { return 0; };

    OneFunction() : Function<dim>(oneFun, zeroFun) {}
 };


#endif //FUNCTION