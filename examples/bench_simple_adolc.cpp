#include <math.h>
#include <cstdlib>
#include <iostream>

#define N 100
#define FILLITER 10

#ifdef ADOLC_TAPELESS
#define NUMBER_DIRECTIONS N
#ifdef OLD_TAPELESS
#include <adouble.h>
#include <adalloc.h>
typedef adtl::adouble adouble;
ADOLC_TAPELESS_UNIQUE_INTERNALS;
#else
#include <adolc/adtl.h>
#include <adalloc.h>
typedef adtl::adouble adouble;
#endif
#else
#include <adolc.h>
#endif

void laplacian(adouble x[N], adouble y[N])
{
    int i, p;
    y[0] = -2*x[0] + x[1];
    y[N-1] = -2*x[N-1] + x[N-2];
    for (i = 1; i < N-1; ++i) {
        y[i] = x[i-1] - 2*x[i] + x[i+1];
    }
    for (p = 0; p < FILLITER; ++p) {
        for (i = 1; i < N-1; ++i) {
            y[i] = y[i-1] - 2*y[i] + y[i+1]/(3.1415+y[i]);
        }
    }
}

#ifdef ADOLC_TAPELESS
void doit(double x[N], double y[N], double **J)
{
#ifndef OLD_TAPELESS
    adtl::setNumDir(N);
#endif

    adouble xad[N], yad[N];
    int i, j;

    for (i = 0; i < N; i++) {
        xad[i] = x[i];
        xad[i].setADValue(i, 1);
    }
    laplacian(xad, yad);
    for (i = 0; i < N; i++) {
        y[i] = yad[i].getValue();
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            J[i][j] = yad[i].getADValue(j);
        }
    }
}
#else
void doit(double x[N], double y[N], double **J)
{
    adouble xad[N], yad[N];
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

    jacobian(1, N, N, x, J);
}
#endif

int main() {
    double x[N];
    double y[N];
    double **J;
    int rep, i, j;

    for (i = 0; i < N; ++i) {
        x[i] = 1.0 + i + 1;
    }

    J = myalloc2(N,N);
    
    for (rep = 0; rep < int(1e8/(N*N)); ++rep) {
        doit(x, y, J);
    }

    for (i = 0; i < 5; ++i) {
        for (j = 0; j < 5; ++j) {
            std::cout << J[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
