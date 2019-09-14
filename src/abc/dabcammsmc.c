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
 * @brief Approximate Bayesian computation moment matching sequential Monte Carlo
 *        with adaptive proposals.
 * @details Likelihood free sequential Monte Carlo using an approximate model
 * to precondition the sampler
 *
 * @param aabc_p ABC Parameters for approximate mode for preconditioning
 * @param abc_p ABC Parameters for exact model
 * @param smc_p SMC parameters
 * @param data dataset to condition on
 * @param atheta array to store resulting particles
 * @param aweights normalised approx  paticle weights
 * @param arho sample discrepency metric values 
 * @param theta array to store resulting exact particles
 * @param weights normalised exact paticle weights
 * @param rho sample discrepency metric values
 *
 * @note current assumption is that the parameterisation of the approximate and 
 * exact model is the same. That is aabc_p and abc_p only differ in the simulation
 * function. Moving forward, it would be interesting to consider extending the 
 * parameter space.
 *
 */
int 
dabcammsmc(ABC_Parameters aabc_p, ABC_Parameters abc_p, SMC_Parameters smc_p,
        Dataset * data, double * atheta, double * aweights, double * arho, 
        double * theta, double * weights,double *rho)
{
    Dataset *data_s;

    double * atheta_prev;
    double * aweights_prev;
    double * arho_prev;
    double aW_sum;
    double aW_sum_prev;
    double * aMu;
    double * aSigma;
    
    double * theta_prev;
    double * weights_prev;
    double * rho_prev;
    double W_sum;
    double W_sum_prev;
    double * Mu;
    double * Sigma;
    clock_t start_t;
    clock_t stamp_t;

    start_t = clock();
    /*allocate memory for simulated  dataset*/
    data_s = copyDataset(data); 
   
    /*allocate for approx samples*/
    atheta_prev = (double *)malloc(aabc_p.k*aabc_p.nacc*sizeof(double));
    aweights_prev = (double *)malloc(aabc_p.nacc*sizeof(double));
    if (arho != NULL)
        arho_prev = (double *)malloc(aabc_p.nacc*sizeof(double));
    aMu = (double *)malloc(aabc_p.k*sizeof(double));
    aSigma = (double *)malloc(aabc_p.k*aabc_p.k*sizeof(double));
    
    /*allocate for exact samples*/
    theta_prev = (double *)malloc(abc_p.k*(abc_p.nacc+aabc_p.nacc)*sizeof(double));
    weights_prev = (double *)malloc((abc_p.nacc+aabc_p.nacc)*sizeof(double));
    if (rho != NULL)
        rho_prev = (double *)malloc((abc_p.nacc+aabc_p.nacc)*sizeof(double));
    Mu = (double *)malloc(abc_p.k*sizeof(double));
    Sigma = (double *)malloc(abc_p.k*abc_p.k*sizeof(double));
    
    /*initialise with ABC rejection samples with eps_0 and approximate model*/
    /*approx*/
    aabc_p.eps = smc_p.eps_t[0];
    dabcrs(aabc_p, data, atheta, arho);
    /*exact pre-sample*/
    abc_p.eps = smc_p.eps_t[0];
    dabcrs(abc_p, data, theta,rho);

   /*initialise weights aW_i = 1 and W_i = 1(we don't normalise, 
    rather we just store the sum)*/
    for (int i=0;i<aabc_p.nacc;i++)
    {
        aweights[i] = 1.0;
    }
    aW_sum = (double)aabc_p.nacc;
    
    for (int i=0;i<abc_p.nacc;i++)
    {
        weights[i] = 1.0;
    }
    W_sum = (double)abc_p.nacc;

    /* compute empirical means and covariances */
    dmcintcv(aabc_p.nacc,aabc_p.k,atheta,aMu,aSigma); /*approximate*/
    dmcintcv(abc_p.nacc,abc_p.k,theta,Mu,Sigma); /*exact*/
    
    /* factorise covariances and transform approx samples*/
    choldc(aabc_p.k,aSigma); /*approximate*/
    choldc(abc_p.k,Sigma);  /*exact*/
    
    /* Apply transform L_exact[(L_approx^-1)(theta_approx - Mu_approx)] + Mu_exact*/
    /* transformed approx samples added to exact pool */
    for (unsigned int i=0;i<aabc_p.nacc;i++)
    {
        double x[abc_p.k];
        for (unsigned int j=0;j<aabc_p.k;j++)
        {
            /*copy approximate samples to the remainder of exact pool*/
            theta[(i+abc_p.nacc)*aabc_p.k + j] = atheta[i*aabc_p.k + j];
            /* de-mean approximate samples*/
            theta[(i+abc_p.nacc)*aabc_p.k + j] -= aMu[j];
        }
        /* solve L_approx (x) = (theta_approx - mu_approx)*/
        cholfs(aabc_p.k, aSigma, x, theta + (i+abc_p.nacc)*aabc_p.k);
        /* compute L_exact (x) */
        mvprod(abc_p.k, Sigma, x, theta + (i+abc_p.nacc)*abc_p.k);
        for (unsigned int j=0;j<aabc_p.k;j++)
        {
            theta[(i+abc_p.nacc)*abc_p.k + j] += Mu[j];
        }
    }
    for (unsigned int i=0;i<aabc_p.nacc;i++)
    {
        weights[i+abc_p.nacc] = aweights[i]; 
    }
    W_sum += aW_sum;
    if (arho != NULL && rho != NULL)
    {
        for (unsigned int i=0;i<aabc_p.nacc;i++)
        {
            rho[i+abc_p.nacc] = arho[i];
        }
    }

#if defined(__CHECKPOINT__)
    {
        /*for long running simulations*/
        /* particle,eps,theta1,...,thetaN,rho?,weights,SEC,N,NSIMS */
        FILE *fp;
        double time;
        fp = fopen(CHECKPOINT_FILENAME,"a");
        stamp_t = clock();
        time = ((double)(stamp_t - start_t))/((double)CLOCKS_PER_SEC);
        /*output particles*/
        for (unsigned int i=0;i<(abc_p.nacc+aabc_p.nacc);i++)
        {
            fprintf(fp,"%d,%lg",i,smc_p.eps_t[0]);
            if (i<aabc_p.nacc)
            {
                /*output pre-transformed particles for diagnositics*/
                for (unsigned int j=0;j<aabc_p.k;j++)
                {
                    fprintf(fp,",%lg",atheta[i*aabc_p.k + j]);
                }
                if (arho != NULL)
                    fprintf(fp,",%lg",arho[i]);
                fprintf(fp,",%lg",aweights[i]);
            }
            else
            {
                for (unsigned int j=0;j<aabc_p.k;j++)
                {
                    fprintf(fp,",0");
                }
                if (arho != NULL)
                    fprintf(fp,",-1");
                fprintf(fp,",-1");
            }
            /*output transformed particles and exact particles*/
            for (unsigned int j=0;j<abc_p.k;j++)
            {
                fprintf(fp,",%lg",theta[i*abc_p.k + j]);
            }
            if (rho != NULL)
                fprintf(fp,",%lg",rho[i]);
            fprintf(fp,",%lg,%lg,%d,%d,%d\n",weights[i],time,abc_p.nacc+aabc_p.nacc,approx_sim_counter,sim_counter);
        }
        fclose(fp);
    }
#endif 
 
    /*commence sequential Monte Carlo steps*/
    for (unsigned int t=1;t<smc_p.T;t++)
    {
        /*copy previous state*/
        
        /*approximate*/
        memcpy(atheta_prev,atheta,aabc_p.k*aabc_p.nacc*sizeof(double));
        memcpy(aweights_prev,aweights,aabc_p.nacc*sizeof(double));
        if (arho != NULL)
            memcpy(arho_prev,arho,aabc_p.nacc*sizeof(double));
        aW_sum_prev = aW_sum;
        
        /*exact*/
        memcpy(theta_prev,theta,abc_p.k*(abc_p.nacc+aabc_p.nacc)*sizeof(double));
        memcpy(weights_prev,weights,(abc_p.nacc+aabc_p.nacc)*sizeof(double));
        if (rho != NULL)
            memcpy(rho_prev,rho,(abc_p.nacc+aabc_p.nacc)*sizeof(double));
        W_sum_prev = W_sum;
       
        /*compute optimal kernel parameters based on particles from 
         * the previous iteration for approximate model*/
        smc_p.adpt(aabc_p.nacc,aabc_p.k,atheta_prev,aweights_prev,smc_p.q_params);
        
        /* importance sampling for approximate model*/
        for (int i=0;i<aabc_p.nacc;i++)
        {
            double d, back_kern;
            d = INFINITY;
            while (d >= smc_p.eps_t[t])
            {
                int j;
                /*sample a particle by weight from {theta_t-1,W_t-1} */
                j = durngpmfs(aabc_p.nacc,aweights_prev,aW_sum_prev);
                /*perturb particle using transition kernel*/
                (*(smc_p.q_adpt))(aabc_p.k,atheta_prev+j*aabc_p.k,
                                           atheta + i*aabc_p.k,
                                           smc_p.q_params);
                /*simulate data with approximate model*/
                (*(aabc_p.s))(aabc_p.sim,atheta +i*aabc_p.k,data_s);
                /*compute discrepency metric*/
                d = (*(aabc_p.rho))(data,data_s);
            }

            if (arho != NULL)
                arho[i] = d;

            /*update particle weight*/
            back_kern = 0;
            for (int j=0;j<aabc_p.nacc;j++)
            {
                back_kern += aweights_prev[j]*((*(smc_p.qd_adpt))(aabc_p.k,
                                                           atheta_prev+j*aabc_p.k,
                                                           atheta+i*aabc_p.k,
                                                           smc_p.q_params));
            }
            aweights[i] = (*(aabc_p.pd))(aabc_p.k,atheta + i*aabc_p.k)*
                         (aW_sum / back_kern);
        }

        /*update weight sum*/
        aW_sum = 0;
        for (int i=0;i<aabc_p.nacc;i++)
        {
            aW_sum += aweights[i];
        }
        
        /*resample with replacement*/
        memcpy(atheta_prev,atheta,aabc_p.k*aabc_p.nacc*sizeof(double));
        if (arho != NULL)
            memcpy(arho_prev,arho,aabc_p.nacc*sizeof(double));
        for (int i=0;i<aabc_p.nacc;i++)
        {
            int j;
            j = durngpmfs(aabc_p.nacc,aweights,aW_sum);
            memcpy(atheta +i*aabc_p.k, atheta_prev + j*aabc_p.k,
                   aabc_p.k*sizeof(double));
            if (arho != NULL)
                arho[i] = arho[j];
        }
        /*reset weights*/
        for (int i=0;i<aabc_p.nacc;i++)
        {
            aweights[i] = 1.0;
        }
        aW_sum = (double)aabc_p.nacc;

        /*compute optimal kernel parameters based on particles from 
         * the previous iteration for exact model*/
        smc_p.adpt(abc_p.nacc+aabc_p.nacc,abc_p.k,theta_prev,weights_prev,smc_p.q_params);
        /* importance sampling for exact model*/
        for (int i=0;i<abc_p.nacc;i++)
        {
            double d, back_kern;
            d = INFINITY;
            while (d >= smc_p.eps_t[t])
            {
                int j;
                /*sample a particle by weight from {theta_t-1,W_t-1} */
                j = durngpmfs(abc_p.nacc+aabc_p.nacc,weights_prev,W_sum_prev);
                /*perturb particle using transition kernel*/
                (*(smc_p.q_adpt))(abc_p.k,theta_prev+j*abc_p.k,
                                          theta + i*abc_p.k,
                                          smc_p.q_params);
                /*simulate data with approximate model*/
                (*(abc_p.s))(abc_p.sim,theta +i*abc_p.k,data_s);
                /*compute discrepency metric*/
                d = (*(abc_p.rho))(data,data_s);
            }

            if (rho != NULL)
                rho[i] = d;

            /*update particle weight*/
            back_kern = 0;
            for (int j=0;j<(abc_p.nacc+aabc_p.nacc);j++)
            {
                back_kern += weights_prev[j]*((*(smc_p.qd_adpt))(abc_p.k,
                                                           theta_prev+j*abc_p.k,
                                                           theta+i*abc_p.k,
                                                           smc_p.q_params));
            }
            weights[i] = (*(abc_p.pd))(abc_p.k,theta + i*abc_p.k)*
                         (W_sum / back_kern);
        }

        /*update weight sum*/
        W_sum = 0;
        for (int i=0;i<abc_p.nacc;i++)
        {
            W_sum += weights[i];
        }
        
        /*resample with replacement*/
        memcpy(theta_prev,theta,abc_p.k*abc_p.nacc*sizeof(double));
        if (rho != NULL)
            memcpy(rho_prev,rho,abc_p.nacc*sizeof(double));
        for (int i=0;i<abc_p.nacc;i++)
        {
            int j;
            j = durngpmfs(abc_p.nacc,weights,W_sum);
            memcpy(theta +i*abc_p.k, theta_prev + j*abc_p.k,
                   abc_p.k*sizeof(double));
            if (rho != NULL)
                rho[i] = rho[j];
        }
        /*reset weights*/
        for (int i=0;i<abc_p.nacc;i++)
        {
            weights[i] = 1.0;
        }
        W_sum = (double)abc_p.nacc;

        /* compute empirical means and covariances */
        dmcintcv(aabc_p.nacc,aabc_p.k,atheta,aMu,aSigma); /*approximate*/
        dmcintcv(abc_p.nacc,abc_p.k,theta,Mu,Sigma); /*exact*/
        
        /* factorise covariances and transform approx samples*/
        choldc(aabc_p.k,aSigma); /*approximate*/
        choldc(abc_p.k,Sigma);  /*exact*/
        
        /* Apply transform L_exact[(L_approx^-1)(theta_approx - Mu_approx)] + Mu_exact*/
        /* transformed approx samples added to exact pool */
        for (unsigned int i=0;i<aabc_p.nacc;i++)
        {
            double x[abc_p.k];
            for (unsigned int j=0;j<aabc_p.k;j++)
            {
                /*copy approximate samples to the remainder of exact pool*/
                theta[(i+abc_p.nacc)*aabc_p.k + j] = atheta[i*aabc_p.k + j];
                /* de-mean approximate samples*/
                theta[(i+abc_p.nacc)*aabc_p.k + j] -= aMu[j];
            }
            /* solve L_approx (x) = (theta_approx - mu_approx)*/
            cholfs(aabc_p.k, aSigma, x, theta + (i+abc_p.nacc)*aabc_p.k);
            /* compute L_exact (x) */
            mvprod(abc_p.k, Sigma, x, theta + (i+abc_p.nacc)*abc_p.k);
            for (unsigned int j=0;j<aabc_p.k;j++)
            {
                theta[(i+abc_p.nacc)*abc_p.k + j] += Mu[j];
            }
        }

        for (unsigned int i=0;i<aabc_p.nacc;i++)
        {
            weights[i+abc_p.nacc] = aweights[i]; 
        }
        W_sum += aW_sum;
        if (arho != NULL && rho != NULL)
        {
            for (unsigned int i=0;i<aabc_p.nacc;i++)
            {
                rho[i+abc_p.nacc] = arho[i];
            }
        }

#if defined(__CHECKPOINT__)
        {
            /*for long running simulations*/
            /* particle,eps,theta1,...,thetaN,rho?,weights,SEC,N,NSIMS */
            FILE *fp;
            double time;
            fp = fopen(CHECKPOINT_FILENAME,"a");
            stamp_t = clock();
            time = ((double)(stamp_t - start_t))/((double)CLOCKS_PER_SEC);
            /*output particles*/
            for (unsigned int i=0;i<(abc_p.nacc+aabc_p.nacc);i++)
            {
                fprintf(fp,"%d,%lg",i,smc_p.eps_t[t]);
                if (i<aabc_p.nacc)
                {
                    /*output pre-transformed particles for diagnositics*/
                    for (unsigned int j=0;j<aabc_p.k;j++)
                    {
                        fprintf(fp,",%lg",atheta[i*aabc_p.k + j]);
                    }
                    if (arho != NULL)
                        fprintf(fp,",%lg",arho[i]);
                    fprintf(fp,",%lg",aweights[i]);
                }
                else
                {
                    for (unsigned int j=0;j<aabc_p.k;j++)
                    {
                        fprintf(fp,",0");
                    }
                    if (arho != NULL)
                        fprintf(fp,",-1");
                    fprintf(fp,",-1");
                }
                /*output transformed particles and exact particles*/
                for (unsigned int j=0;j<abc_p.k;j++)
                {
                    fprintf(fp,",%lg",theta[i*abc_p.k + j]);
                }
                if (rho != NULL)
                    fprintf(fp,",%lg",rho[i]);
                fprintf(fp,",%lg,%lg,%d,%d,%d\n",weights[i],time,abc_p.nacc+aabc_p.nacc,approx_sim_counter,sim_counter);
            }
            fclose(fp);
        }
#endif 
    }

    /*normalise weights*/
    for (int i=0;i<aabc_p.nacc;i++)
    {
        aweights[i] /= aW_sum;
    }
    for (int i=0;i<abc_p.nacc+aabc_p.nacc;i++)
    {
        weights[i] /= W_sum;
    }
    
    free(atheta_prev);     
    free(aweights_prev); 
    if (arho != NULL)
        free(arho_prev);       
    free(aMu);            
    free(aSigma);         
    free(theta_prev);     
    free(weights_prev);   
    if (rho != NULL)
        free(rho_prev);       
    free(Mu);             
    free(Sigma);          
    return 0;
}
