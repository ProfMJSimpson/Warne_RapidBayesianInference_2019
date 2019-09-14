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
 * Computes matrix vector product Ax = b
 */
int
mvprod(unsigned int N, double * A, double *x, double *b)
{
    for (unsigned int i=0;i<N;i++)
    {
        b[i] = 0;
        for (unsigned int k=0;k<N;k++)
        {
            b[i] += A[i*N + k]*x[k];
        }
    }
    return 0;
}
