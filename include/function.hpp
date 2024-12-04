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

    /**
     * @return this a new Function evaluating to its gradient
     */
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

 class ZeroFunction : public Function {
    public:
    ZeroFunction() : Function(0,0) {}

    inline double value(double x) const override{
        return 0.0;
    }
    inline double grad(double x) const override{
        return 0.0;
    }
 };
 
 class OneFunction : public Function {
    public:
    OneFunction() : Function(0,0) {}

    inline double value(double x) const override{
        return 1.0;
    }
    inline double grad(double x) const override{
        return 1.0;
    }
 };


#endif //FUNCTION