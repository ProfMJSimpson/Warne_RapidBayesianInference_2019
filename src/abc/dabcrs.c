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
 * @brief ABC rejection sampling
 * @param abc_p ABC Parameter structure
 * @param data a dataset D to evaluate likelihood 
 * @param theta array of accepted posterior samples
 * @param rho values of rho for each sample
 */
int 
dabcrs(ABC_Parameters abc_p, Dataset * data, double * theta,double *rho)
{
    unsigned int i,j;
    double d;
    Dataset *data_s;
    /*allocate memory for simulated  dataset*/
    data_s = copyDataset(data); 
    /*ABC rejection mode, fixed compute budget nmax > 0, 
     * else fixed acceptance threshold eps*/
    if (abc_p.nmax == 0)
    {
        j=0;
        while (j < abc_p.nacc)
        {
            /*sample the prior and write to rate parameters*/
            (*(abc_p.p))(abc_p.k,1,abc_p.support,theta + j*abc_p.k);
            /*simulate the process*/
            (*(abc_p.s))(abc_p.sim,theta+j*abc_p.k,data_s);

            /*accept/reject based on rho*/
            d = (*(abc_p.rho))(data,data_s);
            
            if (d <= abc_p.eps)
            {
#if defined(__CHECKPOINT__)
                {
                    FILE *fp;
                    fp = fopen(CHECKPOINT_FILENAME,"a");
                    fprintf(fp,"%d",j);
                    for (i=0;i<abc_p.k;i++)
                    {
                        fprintf(fp,",%lg",theta[j*abc_p.k + i]);
                    }
                    if (rho != NULL)
                    {
                        fprintf(fp,",%lg",rho[j]);
                    }
                    fprintf(fp,",%d\n",sim_counter);
                    fclose(fp);
                }
#endif
                /*copy accepted result to output arrays*/
                if (rho != NULL) 
                    rho[j] = d;
                j++;

            }
        }
    }
    else /*fixed compute budget (set predictive prior samples)*/
    {
        /*sample the prior*/
        (*(abc_p.p))(abc_p.k,abc_p.nmax,abc_p.support,theta);
        for (j=0;j<abc_p.nmax;j++)
        {
            /*simulate the process and evaluate rho*/
            (*(abc_p.s))(abc_p.sim,theta+j*abc_p.k,data_s);
            rho[j] = (*(abc_p.rho))(data,data_s);
#if defined(__CHECKPOINT__)
            if ((j+1)%1000 == 0)
            {
                FILE *fp;
                fp = fopen(CHECKPOINT_FILENAME,"a");
                for (int k=(j-999);k<=j;k++)
                {
                    fprintf(fp,"%d",k);
                    for (i=0;i<abc_p.k;i++)
                    {
                        fprintf(fp,",%lg",theta[k*abc_p.k + i]);
                    }
                    fprintf(fp,",%lg,%d\n",rho[k],sim_counter);
                }
                fclose(fp);
            }
#endif

        }
    }
    return 0;
}




