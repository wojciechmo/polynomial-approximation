#ifndef SVD_H
#define SVD_H

#include <vector>
#include <cassert>
#include <cstring>
#include <cmath>
#include "matrix.h"

class SVD
{
public:

    SVD();

    Matrix *U;
    Matrix *S;
    Matrix *VT;
    Matrix *A_pinv;

    int multiply (Matrix *m1, Matrix *m2, Matrix *y);
    Matrix* pinv_compute(Matrix *A, int rows, int cols);
    double* pseudo_inverse_SVD(double* V,double *S,double *U,int rows,int cols);
    void GivensL(double* S_, const size_t dim[2], size_t m,double a, double b);
    void GivensR(double* S_, const size_t dim[2], size_t m, double a, double b);
    void svd_decomp(const size_t dim[2], double* U_, double* S_, double* V_, double eps);
    void svd_compute(int *M, int *N, double *A, int *LDA,double *S, double *U, int *LDU, double *VT, int *LDVT);
};

#endif // SVD_H
