#ifndef FUNCTION
#define FUNCTION

#include <functional>
#include <iostream>
#include <vector>
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

    static const inline fun zeroFun = [](const Point<dim> &p) -> Point<1>
    { return Point<1>::zero(); };

    static const inline fun oneFun = [](const Point<dim> &p) -> Point<1>
    { return Point<1>::one(); };

private:
    fun function;
    std::vector<fun> gradient;

public:
    // Public constructors: necessary to create functions from lambda/callable

    //explicit Function(fun f) : function(f) {}


    // f(x), grad = 0
    explicit Function(fun f) : function(f), gradient{zeroFun} {}

    // 1D: f(x), f_x
    explicit Function(fun f, fun gx) : function(f), gradient{gx} {
        static_assert(dim >= 1, "Function ctor with 1 derivative requires dim>=1");
    }

    // 2D: f(x,y), (f_x, f_y)
    explicit Function(fun f, fun gx, fun gy) : function(f), gradient{gx, gy} {
        static_assert(dim >= 2, "Function ctor with 2 derivatives requires dim>=2");
    }

    // 3D: f(x,y,z), (f_x, f_y, f_z)
    explicit Function(fun f, fun gx, fun gy, fun gz) : function(f), gradient{gx, gy, gz} {
        static_assert(dim >= 3, "Function ctor with 3 derivatives requires dim>=3");
    }

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
        return gradient[0](p);
    }

    inline double dy_value(Point<dim> p) const {
        return gradient[1](p);
    }

    inline double dz_value(Point<dim> p) const {
        return gradient[2](p);
    }

    inline std::vector<double> getGradValues(Point<dim> p) const {
        std::vector<double> out;
        out.reserve(gradient.size());
        for (const fun &g : gradient)
            out.push_back(g(p));
        return out;
    }

    Function<dim, dim> getGrad() const;
};


#include "function.tpp"

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

#endif // FUNCTION