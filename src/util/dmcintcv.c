/* MCL: Monte Carlo Library.
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

#include "mcl.h"
/**
 * @brief Monte Carlo Integration for a random vector
 * @details Approximates E[X]  and Cov[X_i,X_j] for 0<=i<k,i<=j<k
 *  X is a random vectors in R^k.
 *
 * @param N the number of iid samples.
 * @param k dimensionality of X random vector.
 * @param X an Nxk array of iid samples. 
 * @param E output of the expectation (vector of length k).
 * @param Cov output of the co-variance (k x k matrix). Due to symmetry, only lower
 * triangular part is written to
 */
int 
dmcintcv(unsigned int N,unsigned int k, double *X,double *E, double *Cov)
{
    
    for (unsigned int j=0;j<k;j++)
    {
        E[j] = 0;
    }
    for (unsigned int i=0;i<N;i++)
    {
        for (unsigned j=0;j<k;j++)
        {
            E[j] += X[i*k +j];
        }
    }
   
    for (unsigned int j=0;j<k;j++)
    {
        E[j] /= (double)N;
    }

    if (Cov == NULL)
    {
        return 0;
    }

    /* compute variance Cov[X_i,X_j] = E[(X_i - E[X_i])(X_j - E[X_j])]*/
    for (unsigned int j=0;j<k*k;j++)
    {
        Cov[j] = 0;
    }

    for (unsigned int i=0;i<N;i++)
    {
        for (unsigned int j=0;j<k;j++)
        {
            for (unsigned int jj=j;jj<k;jj++)
            {
                Cov[jj*k + j] += X[i*k + j]*X[i*k + jj];
            }
        }
    }

    /* using Bessel's bias correction*/
    for (unsigned int j=0;j<k;j++)
    {
        for (unsigned int jj=j;jj<k;jj++)
        {
            Cov[jj*k + j] /= (double)(N-1);
            Cov[jj*k + j] -= (((double)N)/((double)(N-1)))*E[j]*E[jj];
        }
    }
    return 0;
}
