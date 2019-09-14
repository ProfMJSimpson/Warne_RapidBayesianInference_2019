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
 * Solves the linear system Lx = b where L is a lower triangular matrix
 *
 * @note this function does NOT solve the linear system Ax = b where A is symmetric
 * positive definite. 
 */
int
cholfs(unsigned int N, double * L, double *x, double * b)
{
    for (unsigned int i=0;i<N;i++)
    {
        double sum;
        sum = b[i];
        for (unsigned int k=0;k<i;k++)
        {
            sum -= L[i*N + k]*x[k];
        }
        x[i] = sum/L[i*N + i];
    }
    return 0;
}

