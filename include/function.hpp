#ifndef FUNCTION
#define FUNCTION

#include <functional>

using fun = std::function<double(double)>;

class Function{
    private:
    const fun function;
    const fun gradient;

    public:
    Function(fun f, fun g) : function(f) , gradient(g) {}
    Function(fun f) : function(f) , gradient(0) {}

    inline double value(double x) const{
        return function(x);
    }
    inline double grad(double x) const{
        return gradient(x);
    }

    Function getGrad() const {
        return gradient;
    }

    Function operator +(const Function& f) const;
    Function operator *(const Function& f) const;
    inline double operator ()(double x) const{
        return value(x);
    }

 };


#endif //FUNCTION