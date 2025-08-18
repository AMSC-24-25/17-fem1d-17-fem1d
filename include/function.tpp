#pragma once


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

