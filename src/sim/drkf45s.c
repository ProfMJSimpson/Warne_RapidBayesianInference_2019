/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2019  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "SSAL.h"

/**
 * @brief 4th order Runge-Kutta-Fehlberg (RKF45) method with 5th order error 
 *     estimate and adaptive timesteps.
 * @detail An integration scheme for ordinary differential equations of the form
 *     dY/dt = f(t,Y) where Y \in R^n and f : R^n x [0, \infty) -> R^n
 * @param m     the number of parameters,
 * @param n     dimension of Y
 * @param nt    number of times steps
 * @param T     array of nt time points
 * @param p     vector of model parameters
 * @param Y0    initial condition
 * @param f     vector valued function pointer defining the RHS of the ODE
 * @param ndims number of dimensions to measure
 * @param dims  the dimension indices
 * @param h0    the initial time step
 * @param tol   tolerance on local truncation error
 * @param Y_r   solution stajectory for measured dims (nt*dims)
 */
int 
drkf45s(int m, int n, int nt, double * restrict T, double * restrict p, 
      double * restrict Y0,
      void (*f)(double *, unsigned int, double *,unsigned int, double, double*), 
      int ndims, int * restrict dims, double h0, double tol, 
      double * restrict Y_r)
{
    double Y[n]; /*current state*/
    double Y_i[n]; /*intermediate states*/
    double w[n]; /*update step  */
    double K1[n], K2[n],K3[n], K4[n], K5[n], K6[n]; /*RKF45 intermediate variables*/
    double t = 0; 
    double h = h0;
    double s; /*step scale factor*/
    double epsilon = 2.0*tol;
    int i, ti;
    unsigned int accept = 0;

    for (i=0;i<n;i++)
    {
        Y[i] = Y0[i];
    }
    ti=0;
    while (ti<nt)
    {
        while (!accept)
        {
            double tj;
            double epsilonj;
            /* step 1: K1 = f(t,y)*/
            tj = t;
            (*f)(Y,n,p,m,tj,K1);
            /* step 2: K2 = f(t + h/4,y + h K1/4)*/
            tj = t + h/4.0;
            for(i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + h*K1[i]/4.0;
            }
            (*f)(Y_i,n,p,m,tj,K2);
            /* step 3: K3 = f(t + 3h/8,y +h( 3 K1 /32 + 9 K2 /32))*/
            tj = t + 3.0*h/8.0;
            for (i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + h*(3.0*K1[i] + 9.0*K2[i])/32.0;
            }
            (*f)(Y_i,n,p,m,tj,K3);
            /* step 4: K4 = f(t + 12h/13,y + h(1932 K1 /2197 - 7200 K2 /2197 + 7296 K3 /2197))*/
            tj =t +  12.0*h/13.0;
            for (i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + h*(1932.0*K1[i] - 7200.0*K2[i] + 7296.0*K3[i])/2197.0;
            }
            (*f)(Y_i,n,p,m,tj,K4);
            /* step 5: K5 = f(t + h,y + h(439 K1 /216 - 8 K2 +3680 K3 /513 - 845 K4 /4104)) */
            tj = t + h;
            for (i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + h*(439.0*K1[i]/216.0 -8.0*K2[i] + 3680.0*K3[i]/513.0 - 845.0*K4[i]/2404.0);
            }
            (*f)(Y_i,n,p,m,tj,K5);
            /* step 6: K6 = f(t + h/2, y - h(8 K1 /27 + 2 K2 - 3544 K3 /2565 + 1859 K4 /4104 - 11 K5 /40)) */
            tj = t + h/2.0;
            for(i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + h*(-8.0*K1[i]/27.0 + 2.0*K2[i] - 3544.0*K3[i]/2565.0 + 1859.0*K4[i]/4104.0 - 11.0*K5[i]/40.0);
            }
            (*f)(Y_i,n,p,m,tj,K6);

            /*approximate truncation error*/
            epsilon = fabs(K1[0]/360.0 - 128.0*K3[0]/4275.0 - 2197.0*K4[0]/75240.0 + K5[0]/50.0 + 2.0*K6[0]/55.0);
            for (i=1;i<n;i++)
            {
                epsilonj = fabs(K1[i]/360.0 - 128.0*K3[i]/4275.0 - 2197.0*K4[i]/75240.0 + K5[i]/50.0 + 2.0*K6[i]/55.0);
                epsilon = (epsilonj > epsilon) ? epsilonj : epsilon;
                    
            }

            if (epsilon <= tol)
            {
                /* solution is now within tolerance*/
                accept = 1;
            }
            else /* reject solution, attempt with reduced step */
            {
                
                s = pow(tol/(2.0*epsilon),1.0/4.0);
                s = (s < 0.1) ? 0.1 : ((s > 4.0) ? 4.0 : s);
                h = s*h;
            }
        }
        /*Y update */
        for (i=0;i<n;i++)
        {
            w[i] = (16.0*K1[i]/135.0 + 6656.0*K3[i]/12825.0 + 28561.0*K4[i]/56430.0 - 9.0*K5[i]/50.0 + 2.0*K6[i]/55.0);
        }
        /* check if t < T[ti] < t + h*/
        while (ti < nt && T[ti] < t + h)
        {
            /*interpolate solution and write out timestep*/
            for (i=0;i<ndims;i++)
            {
                
                Y_r[i*nt+ti] = Y[dims[i]] + (T[ti] - t)*(w[dims[i]]);
            }
            ti++;
        }
        /* complete update of solution Y using 5th order extrapolation*/
        for (i=0;i<n;i++)
        {
            Y[i] += h*w[i];
        }
        /* update t*/
        t += h;
        /* update h for next step*/
        s = pow(tol/(2.0*epsilon),1.0/4.0);
        s = (s < 0.1) ? 0.1 : ((s > 4.0) ? 4.0 : s);
        h = s*h;
        accept = 0;
    }
    return 0;
}
