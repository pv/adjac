/*  Example code adapted from adept-1.0 advection example

    Copyright (C) 2012-2013 Robin Hogan and the University of Reading

    Contact email address: r.j.hogan@reading.ac.uk

    This file is part of the Adept library.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <cmath>
#include <iostream>
#include <boost/format.hpp>
#include <vector>

#include <adept.h>

#define NX 100
using namespace adept;

// Lax-Wendroff scheme applied to linear advection
void
lax_wendroff(int nt, double c, const adouble q_init[NX], adouble q[NX])
{
    preallocate_statements((nt+1)*NX*3);
    preallocate_operations((nt+1)*NX*7);
    adouble flux[NX-1];                        // Fluxes between boxes
    for (int i=0; i<NX; i++) {
        q[i] = q_init[i]; // Initialize q
    }
    for (int j=0; j<nt; j++)
    {
        for (int i=0; i<NX-1; i++) {
            flux[i] = 0.5*c*(q[i]+q[i+1]+c*(q[i]-q[i+1]));
        }
        for (int i=1; i<NX-1; i++) {
            q[i] += flux[i-1]-flux[i];
        }
        q[0] = q[NX-2];
        q[NX-1] = q[1];          // Treat boundary conditions
    }
}

// Toon advection scheme applied to linear advection
void
toon(int nt, double c, const adouble q_init[NX], adouble q[NX])
{
    preallocate_statements((nt+1)*NX*3);
    preallocate_operations((nt+1)*NX*9);
    adouble flux[NX-1];                        // Fluxes between boxes
    for (int i=0; i<NX; i++) {
        q[i] = q_init[i]; // Initialize q
    }
    for (int j=0; j<nt; j++) {                 // Main loop in time
        for (int i=0; i<NX-1; i++) {
            flux[i] = (exp(c*log(q[i]/q[i+1]))-1.0) 
                * q[i]*q[i+1] / (q[i]-q[i+1]);
        }
        for (int i=1; i<NX-1; i++) {
            q[i] += flux[i-1]-flux[i];
        }
        q[0] = q[NX-2]; q[NX-1] = q[1];          // Treat boundary
                                                 // conditions
    }
}

int
main(int argc, char** argv)
{
    double pi = 4.0*atan(1.0);
    double q_init[NX];
    
    int nt = 100;
    int nr = 500;
    double dt = 0.125;
    adept::Stack adept_stack;

    for (int i = 0; i < NX; i++) {
        q_init[i] = (0.5+0.5*sin((i*2.0*pi)/(NX-1.5)))+0.0001;
    }

    for (int j = 0; j < nr; j++)
    {
        adouble adept_q_init[NX];
        adouble adept_q[NX];

        adept::set_values(adept_q_init, NX, q_init);
        adept_stack.new_recording();
        toon(nt, dt, adept_q_init, adept_q);

        adept_stack.independent(adept_q_init, NX);
        adept_stack.dependent(adept_q, NX);

        double jacobian[NX*NX];
        adept_stack.jacobian(jacobian);

        if (j == 1) {
            for (int p = 0; p < 5; ++p) {
                std::cout << "  ";
                for (int q = 0; q < 5; ++q) {
                    std::cout << boost::format("%.16e") % jacobian[p + q*NX] << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    return 0;
}
