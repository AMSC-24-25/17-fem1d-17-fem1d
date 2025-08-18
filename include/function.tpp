#pragma once
#ifndef FUNCTION_TPP
#define FUNCTION_TPP

#include "vector.hpp"
#include <stdexcept>

template <unsigned int dim>
Function<dim> Function<dim>::operator+(const Function<dim>& f) const {
    using fun = std::function<double(const Point<dim>&)>;

    const Function<dim>& thisFun = *this;

    fun resultFunction = [thisFun, f](const Point<dim>& p) -> double {
        return thisFun.value(p) + f.value(p);
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1) {
        resultdx = [thisFun, f](const Point<dim>& p) -> double {
            return thisFun.dx_value(p) + f.dx_value(p);
        };
    }
    if constexpr (dim >= 2) {
        resultdy = [thisFun, f](const Point<dim>& p) -> double {
            return thisFun.dy_value(p) + f.dy_value(p);
        };
    }
    if constexpr (dim >= 3) {
        resultdz = [thisFun, f](const Point<dim>& p) -> double {
            return thisFun.dz_value(p) + f.dz_value(p);
        };
    }

    if constexpr (dim == 1) return Function<dim>(resultFunction, resultdx);
    if constexpr (dim == 2) return Function<dim>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3) return Function<dim>(resultFunction, resultdx, resultdy, resultdz); 
}



template<unsigned int dim>
Function<dim> Function<dim>::operator*(const Function<dim> &f) const {
    using fun = std::function<double(const Point<dim>&)>;

    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, f](Point<dim> p) -> double {
        return thisFun.value(p) * f.value(p);
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1) {
        resultdx = [thisFun, f](Point<dim> p) -> double {
            return thisFun.value(p) * f.dx_value(p) + thisFun.dx_value(p) * f.value(p);
        };
    }

    if constexpr (dim >= 2) {
        resultdy = [thisFun, f](Point<dim> p) -> double {
            return thisFun.value(p) * f.dy_value(p) + thisFun.dy_value(p) * f.value(p);
        };
    }

    if constexpr (dim >= 3) {
        resultdz = [thisFun, f](Point<dim> p) -> double {
            return thisFun.value(p) * f.dz_value(p) + thisFun.dz_value(p) * f.value(p);
        };
    }


    if constexpr (dim == 1) return Function<dim>(resultFunction, resultdx);
    if constexpr (dim == 2) return Function<dim>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3) return Function<dim>(resultFunction, resultdx, resultdy, resultdz); 
}

template<unsigned int dim>
Function<dim> Function<dim>::operator+(const double k) const {
    using fun = std::function<double(const Point<dim>&)>;
    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, k](Point<dim> p) -> double {
        return k + thisFun.value(p);
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1) {
        resultdx = [thisFun, k](Point<dim> p) -> double {
            return thisFun.dx_value(p);
        };
    }

    if constexpr (dim >= 2) {
        resultdy = [thisFun, k](Point<dim> p) -> double {
            return thisFun.dy_value(p);
        };
    }

    if constexpr (dim >= 3) {
        resultdz = [thisFun, k](Point<dim> p) -> double {
            return thisFun.dz_value(p);
        };
    }

    if constexpr (dim == 1) return Function<dim>(resultFunction, resultdx);
    if constexpr (dim == 2) return Function<dim>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3) return Function<dim>(resultFunction, resultdx, resultdy, resultdz);

}

template<unsigned int dim>
Function<dim> Function<dim>::operator*(const double k) const {
    using fun = std::function<double(const Point<dim>&)>;
    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, k](Point<dim> p) -> double {
        return k * thisFun.value(p);
    };

    fun resultdx, resultdy, resultdz;

    if constexpr (dim >= 1) {
        resultdx = [thisFun, k](Point<dim> p) -> double {
            return k * thisFun.dx_value(p);
        };
    }

    if constexpr (dim >= 2) {
        resultdy = [thisFun, k](Point<dim> p) -> double {
            return k * thisFun.dy_value(p);
        };
    }

    if constexpr (dim >= 3) {
        resultdz = [thisFun, k](Point<dim> p) -> double {
            return k * thisFun.dz_value(p);
        };
    }

    if constexpr (dim == 1) return Function<dim>(resultFunction, resultdx);
    if constexpr (dim == 2) return Function<dim>(resultFunction, resultdx, resultdy);
    if constexpr (dim == 3) return Function<dim>(resultFunction, resultdx, resultdy, resultdz); 

}

// In-place additions
template <unsigned int dim>
Function<dim>& Function<dim>::operator+=(const Function<dim>& f) {
    *this = (*this + f);
    return *this;
}

template <unsigned int dim>
Function<dim>& Function<dim>::operator+=(double k) {
    *this = (*this + k);
    return *this;
}

// Build a Vector<dim> from gradient component functions
template <unsigned int dim>
Vector<dim> Function<dim>::getGrad() const {
    std::vector<Function<dim>> comps;
    comps.reserve(gradient.size());
    for (const auto& g : gradient) {
        comps.emplace_back(Function(g));
    }
    return Vector<dim>(std::move(comps));
}

// Dot product between two Vector<dim> yielding a scalar Function<dim>
template <unsigned int dim>
Function<dim> dot_product(const Vector<dim>& a, const Vector<dim>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of the same size for dot_product");
    }
    // Accumulate sum of elementwise products
    auto acc = a[0] * b[0];
    for (size_t i = 1; i < a.size(); ++i) {
        acc = acc + (a[i] * b[i]);
    }
    return acc;
}

#endif // FUNCTION_TPP