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
 * @brief Approximate Bayesian computation preconditioned sequential Monte Carlo
 *  with adaptive proposals.
 * @details Likelihood free sequential Monte Carlo using an approximate model
 * to precondition the sampler
 *
 * @param aabc_p ABC Parameters for approximate mode for preconditioning
 * @param abc_p ABC Parameters for exact model
 * @param smc_p SMC parameters
 * @param t_crit critical point to stop precoditioner and begin exact SMC
 * @param data dataset to condition on
 * @param theta array to store resulting particles
 * @param weights normalised paticle weights
 * @param rho sample discrepency metric values
 *
 * @note current assumption is that the parameterisation of the approximate and 
 * exact model is the same. That is aabc_p and abc_p only differ in the simulation
 * function. Moving forward, it would be interesting to consider extending the 
 * parameter space.
 */
int 
dabcapcsmc(ABC_Parameters aabc_p, ABC_Parameters abc_p, SMC_Parameters smc_p, int t_crit,
        Dataset * data, double * theta, double * weights, double *rho)
{
    unsigned int t,i,j;
    Dataset *data_s;
    double *theta_prev;
    double * weights_prev;
    double * rho_prev;
    double W_sum, W_sum_prev;
    clock_t start_t;
    clock_t stamp_t;
    
    t =0;
    start_t = clock();
    
    /*allocate memory for simulated  dataset*/
    data_s = copyDataset(data); 
    
    theta_prev = (double *)malloc(abc_p.k*abc_p.nacc*sizeof(double));
    weights_prev = (double *)malloc(abc_p.nacc*sizeof(double));
    rho_prev = (double *)malloc(abc_p.nacc*sizeof(double));

    /*initialise with ABC rejection samples with eps_0 and approximate model*/
    aabc_p.eps = smc_p.eps_t[0];
    dabcrs(aabc_p,data, theta, rho);
   
   /*initialise weights W_i = 1 (we don't normalise, 
    rather we just store the sum)*/
    for (i=0;i<abc_p.nacc;i++)
    {
        weights[i] = 1.0;
    }
    W_sum = (double)abc_p.nacc;

#if defined(__CHECKPOINT__)
    {
        /*for long running simulations*/
        /* particle,eps,theta1,...,thetaN,rho?,weights,SEC,N,NASIMS,NSIMS */
        FILE *fp;
        double time;
        fp = fopen(CHECKPOINT_FILENAME,"a");
        stamp_t = clock();
        time = ((double)(stamp_t - start_t))/((double)CLOCKS_PER_SEC);
        /*output particles*/
        for (i=0;i<abc_p.nacc;i++)
        {
            fprintf(fp,"%d,%lg",i,smc_p.eps_t[t]);
            for (j=0;j<abc_p.k;j++)
            {
                fprintf(fp,",%lg",theta[i*abc_p.k + j]);
            }
            if (rho != NULL)
                fprintf(fp,",%lg",rho[i]);
            fprintf(fp,",%lg,%lg,%d,%d,%d\n",weights[i],time,abc_p.nacc,approx_sim_counter,sim_counter);
        }
        fclose(fp);
    }
#endif 
 
    /*commence sequential Monte Carlo steps*/
    for (t=1;t<smc_p.T;t++)
    {
        double ESS;
        memcpy(theta_prev,theta,abc_p.k*abc_p.nacc*sizeof(double));
        memcpy(weights_prev,weights,abc_p.nacc*sizeof(double));
        if (rho != NULL)
            memcpy(rho_prev,rho,abc_p.nacc*sizeof(double));
        W_sum_prev = W_sum;

        /*compute optimal kernel parameters based on particles from 
         * the previous iteration*/
        smc_p.adpt(abc_p.nacc,abc_p.k,theta_prev,weights_prev,smc_p.q_params);
        
        /* preconditioning distribution with approxmiate model*/
        for (i=0;i<aabc_p.nacc;i++)
        {
            double d, back_kern;
            d = INFINITY;
            while (d >= smc_p.eps_t[t])
            {
                /*sample a particle by weight from {theta_t-1,W_t-1} */
                j = durngpmfs(aabc_p.nacc,weights_prev,W_sum_prev);
                /*perturb particle using transition kernel*/
                (*(smc_p.q_adpt))(aabc_p.k,theta_prev+j*aabc_p.k,
                                           theta + i*aabc_p.k,smc_p.q_params);
                /*simulate data with approximate model*/
                (*(aabc_p.s))(aabc_p.sim,theta +i*aabc_p.k,data_s);
                /*compute discrepency metric*/
                d = (*(aabc_p.rho))(data,data_s);
            }

            if (rho != NULL)
                rho[i] = d;

            /*update particle weight*/
            back_kern = 0;
            for (j=0;j<aabc_p.nacc;j++)
            {
                back_kern += weights_prev[j]*((*(smc_p.qd_adpt))(aabc_p.k,
                                                           theta_prev+j*aabc_p.k,
                                                           theta+i*aabc_p.k,
                                                           smc_p.q_params));
            }
            weights[i] = (*(aabc_p.pd))(aabc_p.k,theta + i*aabc_p.k)*
                         (W_sum / back_kern);
        }
        /*update weight sum*/
        W_sum = 0;

        for (i=0;i<aabc_p.nacc;i++)
        {
            W_sum += weights[i];
        }
        /*compute effective sample size*/
        ESS = 0;
        for (i=0;i<aabc_p.nacc;i++)
        {
            ESS += weights[i]*weights[i];
        }
        ESS = W_sum*W_sum/ESS;
        if (ESS < smc_p.E) /*resample with replacement*/
        {
            memcpy(theta_prev,theta,aabc_p.k*aabc_p.nacc*sizeof(double));
            if (rho != NULL)
                memcpy(rho_prev,rho,aabc_p.nacc*sizeof(double));
            for (i=0;i<aabc_p.nacc;i++)
            {
                j = durngpmfs(aabc_p.nacc,weights,W_sum);
                memcpy(theta +i*aabc_p.k, theta_prev + j*aabc_p.k,
                       aabc_p.k*sizeof(double));
                if (rho != NULL)
                    rho[i] = rho[j];
            }
            /*reset weights*/
            for (i=0;i<aabc_p.nacc;i++)
            {
                weights[i] = 1.0;
            }
            W_sum = (double)aabc_p.nacc;
        }

        /* importance resampling with exact model*/
        for (i=0;i<abc_p.nacc;i++)
        {
            double d, back_kern;
            d = INFINITY;
            while (d >= smc_p.eps_t[t])
            {
                /*sample a particle by weight from {theta_t-1,W_t-1} */
                j = durngpmfs(abc_p.nacc,weights_prev,W_sum_prev);
                /*perturb particle using transition kernel*/
                (*(smc_p.q))(abc_p.k,theta_prev+j*abc_p.k,
                                          theta + i*abc_p.k);
                /*simulate data with exact model*/
                (*(abc_p.s))(abc_p.sim,theta +i*abc_p.k,data_s);
                /*compute discrepency metric*/
                d = (*(abc_p.rho))(data,data_s);
            }

            if (rho != NULL)
                rho[i] = d;

            /*update particle weight*/
            back_kern = 0;
            for (j=0;j<abc_p.nacc;j++)
            {
                back_kern += weights_prev[j]*((*(smc_p.qd))(abc_p.k,
                                                           theta_prev+j*abc_p.k,
                                                           theta+i*abc_p.k));
            }
            weights[i] = (*(abc_p.pd))(abc_p.k,theta + i*abc_p.k)*
                         (W_sum / back_kern);
        }
        /*update weight sum*/
        W_sum = 0;

        for (i=0;i<abc_p.nacc;i++)
        {
            W_sum += weights[i];
        }
        /*compute effective sample size*/
        ESS = 0;
        for (i=0;i<abc_p.nacc;i++)
        {
            ESS += weights[i]*weights[i];
        }
        ESS = W_sum*W_sum/ESS;
        if (ESS < smc_p.E) /*resample with replacement*/
        {
            memcpy(theta_prev,theta,abc_p.k*abc_p.nacc*sizeof(double));
            if (rho != NULL)
                memcpy(rho_prev,rho,abc_p.nacc*sizeof(double));
            for (i=0;i<abc_p.nacc;i++)
            {
                j = durngpmfs(abc_p.nacc,weights,W_sum);
                memcpy(theta +i*abc_p.k, theta_prev + j*abc_p.k,
                       abc_p.k*sizeof(double));
                if (rho != NULL)
                    rho[i] = rho[j];
            }
            /*reset weights*/
            for (i=0;i<abc_p.nacc;i++)
            {
                weights[i] = 1.0;
            }
            W_sum = (double)abc_p.nacc;
        }

#if defined(__CHECKPOINT__)
        {
            /*for long running simulations*/
            /* particle,eps,theta1,...,thetaN,rho?,weights,SEC,N,NASIMS,NSIMS */
            FILE *fp;
            double time;
            fp = fopen(CHECKPOINT_FILENAME,"a");
            stamp_t = clock();
            time = ((double)(stamp_t - start_t))/((double)CLOCKS_PER_SEC);
            /*output particles*/
            for (i=0;i<abc_p.nacc;i++)
            {
                fprintf(fp,"%d,%lg",i,smc_p.eps_t[t]);
                for (j=0;j<abc_p.k;j++)
                {
                    fprintf(fp,",%lg",theta[i*abc_p.k + j]);
                }
                if (rho != NULL)
                    fprintf(fp,",%lg",rho[i]);
                fprintf(fp,",%lg,%lg,%d,%d,%d\n",weights[i],time,abc_p.nacc,approx_sim_counter,sim_counter);
            }
            fclose(fp);
        }
#endif 
    }

    /*normalise weights*/
    for (i=0;i<abc_p.nacc;i++)
    {
        weights[i] /= W_sum;
    }
    free(theta_prev);
    free(weights_prev);
    if (rho != NULL)
        free(rho_prev);
    return 0;
}
