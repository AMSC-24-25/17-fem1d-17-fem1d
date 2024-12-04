#include "../include/function.hpp"

Function Function::operator +(const Function& f) const{
    const Function& thisFun = *this;

    fun resultFunction = [thisFun, f](double x) -> double {
        return thisFun.value(x) + f.value(x);
    };

    fun resultGradient = [thisFun, f](double x) -> double {
        return thisFun.grad(x) + f.grad(x);
    };
    
    return Function(resultFunction, resultGradient);
}

Function Function::operator *(const Function& f) const{
    const Function& thisFun = *this;

    fun resultFunction = [thisFun, f](double x) -> double {
        return thisFun.value(x) * f.value(x);
    };

    //grad(a*b) = a*grad(b) + grad(a)*b
    fun resultGradient = [thisFun, f](double x) -> double {
        return thisFun.value(x)*f.grad(x) + thisFun.grad(x)*f.value(x);
    };
    
    return Function(resultFunction, resultGradient);
}

Function Function::operator *(const double k) const {
    const Function& thisFun = *this;

    fun resultFunction = [thisFun, k](double x) -> double {
        return k*thisFun.value(x);
    };  

    fun resultGradient = [thisFun, k](double x) -> double {
        return k*thisFun.grad(x);
    };

    return Function(resultFunction, resultGradient);
}
