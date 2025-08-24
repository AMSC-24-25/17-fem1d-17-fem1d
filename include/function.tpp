#ifndef FUNCTION_TPP
#define FUNCTION_TPP
#include "function.hpp"

//----------------general------------

template <unsigned int dim, unsigned int returnDim>
Function<dim, returnDim> Function<dim, returnDim>::operator+(const Function<dim, returnDim> &f) const
{
    using fun = std::function<Point<returnDim>(const Point<dim> &)>;

    const Function<dim, returnDim> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> Point<returnDim>
    {
        return thisFun.value(p) + f.value(p);
    };

    return Function<dim, returnDim>(resultFunction);
}

template <unsigned int dim, unsigned int returnDim>
Function<dim, 1> Function<dim, returnDim>::operator*(const Function<dim, returnDim> &f) const
{
    using fun = std::function<Point<1>(const Point<dim> &)>;

    const Function<dim, returnDim> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> Point<1>
    {
        return thisFun.value(p) * f.value(p);
    };

    return Function<dim, 1>(resultFunction);
}

template <unsigned int dim, unsigned int returnDim>
Function<dim, returnDim> Function<dim, returnDim>::operator+(double val) const
{
    using fun = std::function<Point<returnDim>(const Point<dim> &)>;
    const Function<dim, returnDim> &thisFun = *this;
    Point<returnDim> k(std::vector<double>(returnDim, val));

    fun resultFunction = [thisFun, k](const Point<dim> &p) -> Point<returnDim>
    {
        return thisFun.value(p) + k;
    };

    return Function<dim, returnDim>(resultFunction);
}

template <unsigned int dim, unsigned int returnDim>
Function<dim, returnDim> Function<dim, returnDim>::operator*(double k) const
{
    using fun = std::function<Point<returnDim>(const Point<dim> &)>;
    const Function<dim, returnDim> &thisFun = *this;

    fun resultFunction = [thisFun, k](const Point<dim> &p) -> Point<returnDim>
    {
        return thisFun.value(p) * k;
    };
    return Function<dim, returnDim>(resultFunction);
}

// In-place additions for general template
template <unsigned int dim, unsigned int returnDim>
Function<dim, returnDim> &Function<dim, returnDim>::operator+=(const Function<dim, returnDim> &f)
{
    *this = (*this + f);
    return *this;
}

template <unsigned int dim, unsigned int returnDim>
Function<dim, returnDim> &Function<dim, returnDim>::operator+=(double k)
{
    *this = (*this + k);
    return *this;
}

// ------------------- Specialization for returnDim = 1 --------------------------

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator+(const Function<dim, 1> &f) const
{
    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> double
    {
        return thisFun.value(p) + f.value(p);
    };

    fun_dd resultGrad = [thisFun, f](const Point<dim> &p) -> Point<dim>
    {
        return thisFun.getGradValues(p) + f.getGradValues(p);
    };

    return Function<dim, 1>(resultFunction, resultGrad);
}

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator*(const Function<dim, 1> &f) const
{
    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> double
    {
        return thisFun.value(p) * f.value(p);
    };

    fun_dd resultGrad = [thisFun, f](const Point<dim> &p) -> Point<dim>
    {
        return thisFun.getGradValues(p) * f.value(p) + f.getGradValues(p) * thisFun.value(p);
    };

    return Function<dim, 1>(resultFunction, resultGrad);
}

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator+(double k) const
{
    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, k](const Point<dim> &p) -> double
    {
        return thisFun.value(p) + k;
    };

    return Function<dim, 1>(resultFunction, gradient);
}

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator*(double k) const
{
    using fun = std::function<double(const Point<dim> &)>;
    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, k](const Point<dim> &p) -> double
    {
        return thisFun.value(p) * k;
    };

    fun_dd gradient = [thisFun, k](const Point<dim> &p) -> Point<dim>
    {
        return thisFun.getGradValues(p) * k;
    };

    return Function<dim, 1>(resultFunction, gradient);
}

// In-place additions for Function<dim, 1> specialization
template <unsigned int dim>
Function<dim, 1> &Function<dim, 1>::operator+=(const Function<dim, 1> &f)
{
    *this = (*this + f);
    return *this;
}

template <unsigned int dim>
Function<dim, 1> &Function<dim, 1>::operator+=(double k)
{
    *this = (*this + k);
    return *this;
}

#endif // FUNCTION_TPP