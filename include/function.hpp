#ifndef FUNCTION
#define FUNCTION

#include <functional>
#include <iostream>

using fun = std::function<double(double)>;

const fun zeroFun = [](double x) -> double { return 0; };
const fun oneFun =  [](double x) -> double { return 1; };

class Function{
    private:
    const fun function;
    const fun gradient;

    public:
    explicit Function(fun f, fun g) : function(f) , gradient(g) {}
    explicit Function(fun f) : function(f) , gradient(zeroFun) {}

    inline virtual double value(double x) const{
        return function(x);
    }
    inline virtual double grad(double x) const{
        return gradient(x);
    }

    /**
     * @return this a new Function evaluating to its gradient
     */
    virtual Function getGrad() const {
        return Function(gradient);
    }

    Function operator +(const Function& f) const;
    Function operator *(const Function& f) const;
    Function operator +(const double k) const;
    Function operator *(const double k) const;
    inline virtual double operator ()(double x) const{
        return value(x);
    }

 };

 class ZeroFunction : public Function {
    public:
    ZeroFunction() : Function(zeroFun, zeroFun) {}
 };
 
 class OneFunction : public Function {
    public:
    OneFunction() : Function(oneFun, zeroFun) {}
 };


#endif //FUNCTION