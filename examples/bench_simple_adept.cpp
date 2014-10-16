#include <math.h>
#include <cstdlib>
#include <iostream>

#include <adept.h>

#define N 100
#define FILLITER 10

using adept::adouble;

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
            y[i] = y[i-1] - 2*y[i] + y[i+1]/(3.1415 + y[i]);
        }
    }
}

void doit(double *x, double *y, double *J)
{
    adept::Stack stack;
    adouble xad[N], yad[N];
    int i;

    adept::set_values(xad, N, x);
    stack.new_recording();

    laplacian(xad, yad);

    stack.independent(xad, N);
    stack.dependent(yad, N);
    stack.jacobian(J);
}


int main() {
    double x[N];
    double y[N];
    double J[N*N];
    int rep, i, j;

    for (i = 0; i < N; ++i) {
        x[i] = 1.0 + i + 1;
    }

    for (rep = 0; rep < int(1e8/(N*N)); ++rep) {
        doit(x, y, J);
    }

    for (i = 0; i < 5; ++i) {
        std::cout << y[i] << std::endl;
    }
    for (i = 0; i < 5; ++i) {
        for (j = 0; j < 5; ++j) {
            std::cout << J[i + N*j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
