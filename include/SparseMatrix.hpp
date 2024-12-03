//
// Created by federico/daniele on 03/12/24.
//
#include <vector>

#ifndef INC_17_FEM1D_17_FEM1D_SPARSEMATRIX_H
#define INC_17_FEM1D_17_FEM1D_SPARSEMATRIX_H


using std::vector;

class SparseMatrix {

private:
    std::vector<double> values;
    vector<int> col_index;
    vector<int> beginner_value;
    int nnz;

public:
    SparseMatrix();

    int nnz(){
        return nnz;
    }

    void incrementNnz(){
        nnz ++;
    }

    void addElement();

};


#endif //INC_17_FEM1D_17_FEM1D_SPARSEMATRIX_H
