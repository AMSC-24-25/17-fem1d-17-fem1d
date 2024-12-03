#include "function.hpp"


Function Function::operator +(const Function& f) const{
    fun resultFunction = [this, f](double x) -> double {
        return this->value(x) + f.value(x);
    };
    
    return Function(resultFunction);
}
Function Function::operator *(const Function& f) const{
    fun resultFunction = [this, f](double x) -> double {
        return this->value(x) * f.value(x);
    };
    
    return Function(resultFunction);
}