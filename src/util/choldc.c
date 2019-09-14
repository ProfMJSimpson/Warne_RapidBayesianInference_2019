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
 * Computes Cholesky decomposition of a symmetric positive definite matrix A = L*L'
 * @note assumes only lower triangular part is stored. 
 * Decomposition is stored in-place.
 *
 * @param N dimension of matrix 
 * @param A lower diagonal of positive definite matrix, L is written to A
 */
int
choldc(unsigned int N, double * A)
{
    for (unsigned int i=0;i<N;i++)
    {
        for (unsigned int j=0;j<i;j++)
        {
            for (unsigned int k=0;k<j;k++)
            {
                A[i*N + j] -= A[i*N + k]*A[j*N + k];
            }
            A[i*N + j] /= A[j*N + j];
        }

        for (unsigned int k=0;k<i;k++)
        {
            A[i*N + i] -= A[i*N + k]*A[i*N + k];
        }
        
        if (A[i*N+i] <= 0.0)
        {
            return 1; /*not positive definite!!*/
        }
        A[i*N + i] = sqrt(A[i*N + i]);
    }
    return 0;
}

