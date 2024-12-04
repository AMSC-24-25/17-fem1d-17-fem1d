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

    inline virtual double value(double x) const{
        return function(x);
    }
    inline virtual double grad(double x) const{
        return gradient(x);
    }

    Function getGrad() const {
        return gradient;
    }

    Function operator +(const Function& f) const;
    Function operator *(const Function& f) const;
    Function operator +(const double k) const;
    Function operator *(const double k) const;
    inline double operator ()(double x) const{
        return value(x);
    }

 };


#endif //FUNCTION