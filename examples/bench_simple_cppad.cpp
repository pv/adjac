#include <math.h>
#include <cstdlib>
#include <iostream>
#include <vector>

#define N 100
#define FILLITER 10

#include <cppad/cppad.hpp>

using CppAD::AD;

void laplacian(AD<double> x[N], AD<double> y[N])
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

void doit(std::vector<double> &x, std::vector<double> &y, std::vector<double> &J)
{
    std::vector< AD<double> > xad(N), yad(N);
    int i;

    for (i = 0; i < N; i++) {
        xad[i] = x[i];
    }
    CppAD::Independent(xad);

    laplacian(xad.data(), yad.data());

    CppAD::ADFun<double> f(xad, yad);
    J = f.Jacobian(x);
}

int main() {
    std::vector<double> x(N), y(N);
    std::vector<double> J;
    int rep, i, j;

    for (i = 0; i < N; ++i) {
        x[i] = 1.0 + i + 1;
    }

    for (rep = 0; rep < 1 + int(1e8/(N*N)); ++rep) {
        doit(x, y, J);
    }

    for (i = 0; i < 5; ++i) {
        for (j = 0; j < 5; ++j) {
            std::cout << J[i*N + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
