/* MCL: Monte Carlo Library.
 * Copyright (C) 2017  David J. Warne
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

#include "mcl.h"
/**
 * @brief Monte Carlo Integration for a random vector
 * @details Approximates E[X]  and V[X] = E[|X - E[X]|^2_\infty]
 *  X is a random vectors in R^k.
 *
 * @param N the number of iid samples.
 * @param k dimensionality of X random vector.
 * @param X an Nxk array of iid samples. 
 * @param E output of the expectation (vector of length k).
 * @param V output of the variance (scalar).
 */
int 
dmcintv(unsigned int N,unsigned int k, double *X,double *E, double *V)
{
    unsigned int i,j;
    double fX_i,fY_i;
    
    for (j=0;j<k;j++)
    {
        E[j] = 0;
    }
    for (i=0;i<N;i++)
    {
        for (j=0;j<k;j++)
        {
            E[j] += X[i*k +j];
        }
    }
   
    for (j=0;j<k;j++)
    {
        E[j] /= (double)N;
    }

    if (V == NULL)
    {
        return 0;
    }

    /* compute variance Var[X] = E[|X - E[X]|_\infty^2]*/
    for (j=0;j<k;j++)
    {
        *V = 0;
    }

    for (i=0;i<N;i++)
    {
        double norm_inf,d;
        norm_inf = 0;
        for (j=0;j<k;j++)
        {
            d = fabs(X[i*k + j] - E[j]);
            norm_inf = (norm_inf < d) ? d : norm_inf;
        }
        *V += norm_inf;
    }
    *V /= (double)N;
    return 0;
}
