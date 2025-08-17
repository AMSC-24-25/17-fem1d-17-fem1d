#include "function.hpp"


template<unsigned int dim>
Function<dim> Function<dim>::operator+(const Function<dim> &f) const {
    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, f](Point<dim> p) -> double {
        return thisFun.value(p) + f.value(p);
    };

    fun resultdx = [thisFun, f](Point<dim> p) -> double {
        return thisFun.dx(p) + f.dx(p);
    };

    fun resultdy = [thisFun, f](Point<dim> p) -> double {
        return thisFun.dy(p) + f.dy(p);
    };

    fun resultdz = [thisFun, f](Point<dim> p) -> double {
        return thisFun.dz(p) + f.dz(p);
    };

    std::vector<fun> resultGradient = {resultdx, resultdy, resultdz};

    return Function<dim>(resultFunction, resultGradient);
}


template<unsigned int dim>
Function<dim> Function<dim>::operator*(const Function<dim> &f) const {
    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, f](Point<dim> p) -> double {
        return thisFun.value(p) * f.value(p);
    };

    fun resultdx = [thisFun, f](Point<dim> p) -> double {
        return thisFun.value(p) * f.dx_value(p) + thisFun.dx_value(p) * f.value(p);
    };

    fun resultdy = [thisFun, f](Point<dim> p) -> double {
        return thisFun.value(p) * f.dy_value(p) + thisFun.dy_value(p) * f.value(p);
    };

    fun resultdz = [thisFun, f](Point<dim> p) -> double {
        return thisFun.value(p) * f.dz_value(p) + thisFun.dz_value  (p) * f.value(p);
    };

    std::vector<fun> resultGradient = {resultdx, resultdy, resultdz};
    return Function<dim>(resultFunction, resultGradient);
}

template<unsigned int dim>
Function<dim> Function<dim>::operator+(const double k) const {
    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, k](Point<dim> p) -> double {
        return k + thisFun.value(p);
    };

    fun resultdx = [thisFun, k](Point<dim> p) -> double {
        return thisFun.dx_value(p);
    };

    fun resultdy = [thisFun, k](Point<dim> p) -> double {
        return thisFun.dy_value(p);
    };

    fun resultdz = [thisFun, k](Point<dim> p) -> double {
        return thisFun.dz_value(p);
    };

    std::vector<fun> resultGradient = {resultdx, resultdy, resultdz};

    return Function<dim>(resultFunction, resultGradient);
}

template<unsigned int dim>
Function<dim> Function<dim>::operator*(const double k) const {
    const Function<dim> &thisFun = *this;

    fun resultFunction = [thisFun, k](Point<dim> p) -> double {
        return k * thisFun.value(p);
    };

    fun resultdx = [thisFun, k](Point<dim> p) -> double {
        return k * thisFun.dx_value(p);
    };

    fun resultdy = [thisFun, k](Point<dim> p) -> double {
        return k * thisFun.dy_value(p);
    };

    fun resultdz = [thisFun, k](Point<dim> p) -> double {
        return k * thisFun.dz_value(p);
    };

    std::vector<fun> resultGradient = {resultdx, resultdy, resultdz};

    return Function<dim>(resultFunction, resultGradient);
}

