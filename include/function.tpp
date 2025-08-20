#pragma once
#ifndef FUNCTION_TPP
#define FUNCTION_TPP

#include "vector.hpp"
#include <stdexcept>

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator+(const Function<dim, 1> &f) const
{
    using fun = std::function<double(const Point<dim> &)>;
    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> double
    {
        // TODO: somma componente per componente se Point<returnDim> non supporta operator+
        return thisFun.value(p) + f.value(p);
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1)
    {
        resultdx = [thisFun, f](const Point<dim> &p) -> double
        {
            // TODO: usa derivate vettoriali (Jacobian) se disponibili
            return thisFun.dx_value(p) + f.dx_value(p);
        };
    }
    if constexpr (dim >= 2)
    {
        resultdy = [thisFun, f](const Point<dim> &p) -> double
        {
            return thisFun.dy_value(p) + f.dy_value(p);
        };
    }
    if constexpr (dim >= 3)
    {
        resultdz = [thisFun, f](const Point<dim> &p) -> double
        {
            return thisFun.dz_value(p) + f.dz_value(p);
        };
    }

    if constexpr (dim == 1)
        return Function<dim, 1>(resultFunction, resultdx);
    if constexpr (dim == 2)
        return Function<dim, 1>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3)
        return Function<dim, 1>(resultFunction, resultdx, resultdy, resultdz);
}

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator*(const Function<dim, 1> &f) const
{
    using fun = std::function<double(const Point<dim> &)>;

    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> double
    {
        // TODO: prodotto componente per componente se Point<returnDim> non supporta operator* con Point
        return thisFun.value(p) * f.value(p);
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1)
    {
        resultdx = [thisFun, f](const Point<dim> &p) -> double
        {
            return thisFun.value(p) * f.dx_value(p) + thisFun.dx_value(p) * f.value(p);
        };
    }
    if constexpr (dim >= 2)
    {
        resultdy = [thisFun, f](const Point<dim> &p) -> double
        {
            return thisFun.value(p) * f.dy_value(p) + thisFun.dy_value(p) * f.value(p);
        };
    }
    if constexpr (dim >= 3)
    {
        resultdz = [thisFun, f](const Point<dim> &p) -> double
        {
            return thisFun.value(p) * f.dz_value(p) + thisFun.dz_value(p) * f.value(p);
        };
    }

    if constexpr (dim == 1)
        return Function<dim, 1>(resultFunction, resultdx);
    if constexpr (dim == 2)
        return Function<dim, 1>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3)
        return Function<dim, 1>(resultFunction, resultdx, resultdy, resultdz);
}

template <unsigned int dim>
Function<dim, 1> Function<dim, 1>::operator+(double k) const
{
    using fun = std::function<double(const Point<dim> &)>;
    const Function<dim, 1> &thisFun = *this;

    fun resultFunction = [thisFun, k](const Point<dim> &p) -> double
    {
        return thisFun.value(p) + k;
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1)
    {
        resultdx = [thisFun, k](const Point<dim> &p) -> double
        {
            return thisFun.dx_value(p);
        };
    }
    if constexpr (dim >= 2)
    {
        resultdy = [thisFun, k](const Point<dim> &p) -> double
        {
            return thisFun.dy_value(p);
        };
    }
    if constexpr (dim >= 3)
    {
        resultdz = [thisFun, k](const Point<dim> &p) -> double
        {
            return thisFun.dz_value(p);
        };
    }

    if constexpr (dim == 1)
        return Function<dim, 1>(resultFunction, resultdx);
    if constexpr (dim == 2)
        return Function<dim, 1>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3)
        return Function<dim, 1>(resultFunction, resultdx, resultdy, resultdz);
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

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1)
    {
        resultdx = [thisFun, k](const Point<dim> &p) -> double
        {
            return thisFun.dx_value(p) * k;
        };
    }
    if constexpr (dim >= 2)
    {
        resultdy = [thisFun, k](const Point<dim> &p) -> double
        {
            return thisFun.dy_value(p) * k;
        };
    }
    if constexpr (dim >= 3)
    {
        resultdz = [thisFun, k](const Point<dim> &p) -> double
        {
            return thisFun.dz_value(p) * k;
        };
    }

    if constexpr (dim == 1)
        return Function<dim, 1>(resultFunction, resultdx);
    if constexpr (dim == 2)
        return Function<dim, 1>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3)
        return Function<dim, 1>(resultFunction, resultdx, resultdy, resultdz);
}

template <unsigned int dim>
FunctionVector<dim> Function<dim, 1>::getGrad() const
{
    std::vector<Function<dim, 1>> comps;
    comps.reserve(gradient.size());
    for (const auto &g : gradient)
    {
        comps.emplace_back(Function<dim, 1>(g));
    }
    return FunctionVector<dim>(std::move(comps));
}

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
Function<dim, returnDim> Function<dim, returnDim>::operator*(const Function<dim, returnDim> &f) const
{
    using fun = std::function<Point<returnDim>(const Point<dim> &)>;

    const Function<dim, returnDim> &thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim> &p) -> Point<returnDim>
    {
        return thisFun.value(p) * f.value(p);
    };

    return Function<dim, returnDim>(resultFunction);
}

template <unsigned int dim, unsigned int returnDim>
Function<dim, returnDim> Function<dim, returnDim>::operator+(double k) const
{
    using fun = std::function<Point<returnDim>(const Point<dim> &)>;
    const Function<dim, returnDim> &thisFun = *this;

    fun resultFunction = [thisFun, k](const Point<dim> &p) -> Point<returnDim>
    {
        auto v = thisFun.value(p);
        for (unsigned int i = 0; i < returnDim; ++i)
            v[i] += k;
        return v;
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
        auto v = thisFun.value(p);
        for (unsigned int i = 0; i < returnDim; ++i)
            v[i] *= k;
        return v;
    };
    return Function<dim, returnDim>(resultFunction);
}

// In-place additions
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

#endif // FUNCTION_TPP