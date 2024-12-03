#include "../include/grid1D.hpp"


FunctionVector Grid1D::getPhiFunctions() const{
    //phi_i = 1 on x_i ; 0 on x_j (j != i)
    FunctionVector phiFunctions;
    
    fun phi0 = [this](double x) -> double {
        if(x < start) return 0.0; //exception??
        if(x < start+h) return 1 - (x-start)/h;
        if(x >= start+h) return 0.0;
    };

    fun phiN = [this](double x) -> double {
        if (x < end -h) return 0.0;
        if (x < end) return (end-x)/h;
        if (x >= end) return 0.0; //exception??
    };
    
    phiFunctions.reserve(N+1);
    phiFunctions.push_back(phi0);
    for (int k = 1; k < N-1; k++)
    {
        double xk = (*this)(k);
        fun phi = [this, xk](double x) -> double {
            if(x <= xk-h) return 0.0;
            if(x < xk) return (x - xk-h)/h;
            if(x < xk+h) return 1 - (x-xk)/h;
            if(x >= xk+h) return 0.0;
        };

        fun gradPhi = [this, xk](double x) -> double {
            if(x <= xk-h) return 0.0;
            if(x < xk) return N;
            if(x < xk+h) return -N;
            if(x >= xk+h) return 0.0;
        };
        
        phiFunctions.push_back(Function(phi, gradPhi));
    }
    phiFunctions.push_back(phiN);
    
    return phiFunctions; 
}