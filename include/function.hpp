#ifndef FUNCTION
#define FUNCTION

#include <functional>

using fun = std::function<double(double)>;

class Function{
    public:
    
    Function(fun _f) : function(_f) {}
    
    inline double Function::value(double x) const{
        return function(x);
    }
    inline double Function::grad(double x) const{
        return gradient(x);
    }

    Function operator +(const Function& f) const;
    Function operator *(const Function& f) const;
    inline double operator ()(double x) const{
        return value(x);
    }

    private:

    fun function;
    fun gradient;
 };


#endif //FUNCTION