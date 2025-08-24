#ifndef FUNCTION
#define FUNCTION

#include <functional>
#include <iostream>
#include <vector>
#include <stdexcept>
#include "point.hpp"

template<unsigned int dim, unsigned int returnDim>
using fun_td = std::function<Point<returnDim>(const Point<dim>&, double)>;

template <unsigned int dim, unsigned int returnDim>
class Function
{
    using fun = std::function<Point<returnDim>(const Point<dim> &)>;

public:
    static const inline fun zeroFun = [](const Point<dim>&) -> Point<returnDim> {
        return Point<returnDim>::zero();
    };
    static const inline fun oneFun = [](const Point<dim>&) -> Point<returnDim> {
        return Point<returnDim>::one();
    };

private:
    fun function;

public:

    explicit Function(fun f) : function(f) {}

    Function<dim, returnDim> operator+(const Function<dim, returnDim> &f) const;
    Function<dim, 1> operator*(const Function<dim, returnDim> &f) const;
    Function<dim, returnDim> operator+(double k) const;
    Function<dim, returnDim> operator*(double k) const;
    Function<dim, returnDim> &operator+=(const Function<dim, returnDim> &f);
    Function<dim, returnDim> &operator+=(double k);

    Point<returnDim> value(const Point<dim> &p) const {
        return function(p);
    }

    Point<returnDim> operator()(const Point<dim> &p) const {
        return function(p);
    }

    operator fun_td<dim, returnDim>() const {
        return [f = this->function](const Point<dim>& p, double) -> Point<returnDim> {
            return f(p);
        };
    }
};

template<unsigned int dim, unsigned int returnDim>
inline Function<dim, returnDim> castToSteadyFunction(const fun_td<dim, returnDim>& f_td, double t) {
    return Function<dim, returnDim>([f_td, t](const Point<dim>& p) -> Point<returnDim> {
        return f_td(p, t);
    });
}

template<unsigned dim>
class Function<dim, 1>{

    using fun = std::function<Point<1>(const Point<dim> &)>;
    using fun_dd = std::function<Point<dim>(const Point<dim> &)>;

    static const inline fun zeroFun = [](const Point<dim> &p) -> Point<1>
    { return Point<1>::zero(); };

    static const inline fun oneFun = [](const Point<dim> &p) -> Point<1>
    { return Point<1>::one(); };

private:
    fun function;
    fun_dd gradient;

public:
    // Public constructors: necessary to create functions from lambda/callable

    //explicit Function(fun f) : function(f) {}


    // f(x), grad = 0
    explicit Function(fun f) : function(f), gradient{Function<dim,dim>::zeroFun} {}
    explicit Function(fun f, fun_dd g) : function(f), gradient{g} {}

public: 

    Function<dim, 1> operator+(const Function<dim, 1> &f) const;
    Function<dim, 1> operator*(const Function<dim, 1> &f) const;
    Function<dim, 1> operator+(double k) const;
    Function<dim, 1> operator*(double k) const;
    Function<dim, 1> &operator+=(const Function<dim, 1> &f);
    Function<dim, 1> &operator+=(double k);

    double value(Point<dim> p) const {
        return function(p);
    }

    double operator()(Point<dim> p) const {
        return function(p);
    }

    inline double dx_value(Point<dim> p) const {
        return gradient(p)[0];
    }

    inline double dy_value(Point<dim> p) const {
        return gradient(p)[1];
    }

    inline double dz_value(Point<dim> p) const {
        return gradient(p)[2];
    }

    inline Point<dim> getGradValues(Point<dim> p) const {
        return gradient(p);
    }

    inline Function<dim, dim> getGrad() const {
        return Function<dim, dim>(gradient);
    }
};

template <unsigned int dim, unsigned int returnDim>
class ZeroFunction : public Function<dim, returnDim>
{

    using fun = std::function<Point<returnDim>(const Point<dim> &)>;

public:
    static const inline fun zeroFun = [](const Point<dim> &p) -> Point<returnDim>
    { return Point<returnDim>::zero(); };

    ZeroFunction() : Function<dim, returnDim>(zeroFun, zeroFun) {}
};

template <unsigned int dim, unsigned int returnDim>
class OneFunction : public Function<dim, returnDim>
{

    using fun = std::function<Point<returnDim>(const Point<dim> &)>;

public:
    static const inline fun oneFun = [](const Point<dim> &p) -> Point<returnDim>
    { return Point<returnDim>::one(); };
    static const inline fun zeroFun = [](const Point<dim> &p) -> Point<returnDim>
    { return Point<returnDim>::zero(); };

    OneFunction() : Function<dim, returnDim>(oneFun, zeroFun) {}
};

#include "function.tpp"

#endif // FUNCTION