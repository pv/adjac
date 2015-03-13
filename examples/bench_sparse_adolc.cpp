#include <math.h>
#include <cstdlib>
#include <iostream>

#define BLOCKSIZE 128
#define N (1000*BLOCKSIZE)
#define FILLITER 3

#include <adolc.h>
#include <adolc_sparse.h>

void laplacian(adouble x[N], adouble y[N])
{
    int i, p, k;
    y[0] = -2*x[0] + x[1];
    y[N-1] = -2*x[N-1] + x[N-2];
    for (i = 1; i < N-1; ++i) {
        y[i] = x[i-1] - 2*x[i] + x[i+1];
    }
    for (k = 0; k < N; k += BLOCKSIZE) {
        for (p = 0; p < FILLITER; ++p) {
            for (i = k+1; i < k+BLOCKSIZE-1; ++i) {
                y[i] = y[i-1] - 2*y[i] + y[i+1]/(3.1415+y[i]);
            }
        }
    }
}

void doit(double x[N], double y[N], int *nnz, unsigned int **ii, unsigned int **jj, double **vv)
{
    adouble xad[N], yad[N];
    int options[4] = {0,0,0,0};
    int i;

    trace_on(1);
    for (i = 0; i < N; i++) {
        xad[i] <<= x[i];
    }
    laplacian(xad, yad);
    for (i = 0; i < N; i++) {
        yad[i] >>= y[i];
    }
    trace_off();

    options[0] = 0;
    sparse_jac(1, N, N, 0, x, nnz, ii, jj, vv, options);
}

int main() {
    double x[N];
    double y[N];

    int nnz;
    unsigned int *ii = NULL, *jj = NULL;
    double *vv = NULL;
    int rep, i, j;

    for (i = 0; i < N; ++i) {
        x[i] = 1.0 + i + 1;
    }

    for (rep = 0; rep < 5; ++rep) {
        free(ii);
        free(jj);
        free(vv);
        ii = NULL;
        jj = NULL;
        vv = NULL;
        doit(x, y, &nnz, &ii, &jj, &vv);
    }

    for (i = 0; i < 5; ++i) {
        std::cout << ii[i] << " ";
    }
    std::cout << std::endl;
    for (i = 0; i < 5; ++i) {
        std::cout << jj[i] << " ";
    }
    std::cout << std::endl;
    for (i = 0; i < 5; ++i) {
        std::cout << vv[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
