#include "../include/function.hpp"

Function Function::operator +(const Function& f) const{
    fun resultFunction = [this, f](double x) -> double {
        return this->value(x) + f.value(x);
    };

    fun resultGradient = [this, f](double x) -> double {
        return this->grad(x) + f.grad(x);
    };
    
    return Function(resultFunction, resultGradient);
}

Function Function::operator *(const Function& f) const{
    fun resultFunction = [this, f](double x) -> double {
        return this->value(x) * f.value(x);
    };

    //grad(a*b) = a*grad(b) + grad(a)*b
    fun resultGradient = [this, f](double x) -> double {
        return this->value(x)*f.grad(x) + this->grad(x)*f.value(x);
    };
    
    return Function(resultFunction, resultGradient);
}
