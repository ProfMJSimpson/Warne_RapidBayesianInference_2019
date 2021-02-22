/* LIBABC: approximate Bayesian Computation
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
#include <time.h>

/* Demonstration example of ABC Preconditioned SMC inference using the 
 * Orstein-Uhlenbeck SDE model with the stationary distribution
 * for the preconditioner
 */
ABC_Parameters aabc_p; /* for the continuum model*/
ABC_Parameters abc_p; /* for the discrete model*/
SMC_Parameters smc_p;


/* data set structure */
Dataset data;

/* fixed parameters for the OU model*/
double gam;
double mu;
double T;
double dt; 
double x0;

/**
 * @brief Discrepency metric
 * @details Computes the discrepancy measure \rho(D,D_s) used for ABC inference.
 *          For the weak Alleey model the Euclidean distance is used.
 *
 * @param data      pointer to data used for inference
 * @param data_s    pointer to simulated data generated by forward model
 * 
 * @return the descrepancy measure for thresholding in ABC
 */
double 
rho(Dataset *data, Dataset *data_s)
{
    double *m_s,*m_d; /* data process  generates moments */
    double d;
    
    m_s = (double *)data_s->fields[0].data_array;
    m_d = (double *)data->fields[0].data_array;
    /* rho = |mu_d - m_s| + |sig_d - sig_s|*/
    //d = fabs(m_d[0] - m_s[0]) + fabs(sqrt(m_d[1]) - sqrt(m_s[1]));
    d = fabs(sqrt(m_d[1]) - sqrt(m_s[1]));
    return d;
}

/**
 * @brief prior sampler 
 * @details Generates a number of i.i.d. samples from p(D)
 *
 * @param dim       dimensionality of parameter space 
 * @param nsamples  number of i.i.d. samplesto generate
 * @param support   pointer to array of 2*dim elements containing support limits
 * @param theta     pointer to nsamples*dim array to store samples.
 */
int 
prior(unsigned int dim, unsigned int nsamples, double* support,double * theta)
{
    unsigned int i,j;
    for (i=0;i<nsamples;i++)
    {
        theta[i] = durngus(support[0],support[1]); // D 
    }
    return 0;
}

/**
 * @brief prior probability density function
 * @details evaluates the prior density pointwise
 *
 * @param dim   dimensionality of parameter space
 * @param theta pointer to point to evaluate the PDF at
 *
 * @returns the value of the prior density at theta
 */
double
priorPDF(unsigned int dim, double *theta)
{
    unsigned char tf;
    double p_inv;
    tf = 1; 
    p_inv = 1.0;
    for (int j=0;j<dim;j++)
    {
        tf = tf && (theta[j] >= abc_p.support[j]) && (theta[j] <= abc_p.support[j+dim]);
        p_inv *= (abc_p.support[j+dim] - abc_p.support[j]);
    }
    return (tf) ? 1.0/p_inv : 0.0;
}

#define SCALE 2.0

/**
 * @brief Adaptive updates too proposal kernel parameters
 * @details Uses particle population at step t-1 to compute parameters for
 *           proposal kernels use at step t. In this example, for Gaussian 
 *           proposals, we take the proposal covariance to be twice the particle
 *           covariance at t-1.
 *
 * @param N         the number of particles
 * @param dim       dimensionality if parameter space
 * @param theta     particles from SMC step t-1
 * @param weights   particle weights from step t-1
 * @param L         pointer to store kernel parameter values. In this case it is 
 *                  a dim*dim array to store the lower triangular matrix from 
 *                  the Cholesky factorization of the covariance matrix.
 */
int
adapt_kernel(unsigned int N, unsigned int dim, double * theta, double *weights, 
             double *L)
{
    double mu[dim], weight_bar;
    unsigned int i,j,k;
    for (k=0;k<dim;k++)
    {
        mu[k] = 0.0;
    }
    weight_bar = 0;
    for (k=0;k<dim*dim;k++)
    {
        L[k] = 0.0;
    }

    /* compute weighted mean*/
    for (i=0;i<N;i++)
    {
        for (k=0;k<dim;k++)
        {
            mu[k] += weights[i]*theta[i*dim + k];
        }
        weight_bar += weights[i];
    }
    for (k=0;k<dim;k++)
    {
        mu[k] /= weight_bar;
    }

    /* unbiased estimate of sample covariance (store in L though)*/
    for (i=0;i<N;i++)
    {
        for (j=0;j<dim;j++)
        {
            for (k=0;k<dim;k++)
            {
                L[j*dim + k] += weights[i]*(theta[i*dim + j] - mu[j])
                                          *(theta[i*dim + k] - mu[k]); 
            }
        }
    }
    for (k=0;k<dim*dim;k++)
    {
        L[k] /= (weight_bar - 1.0);
    }

    /*optimal scale is 2 based on the work of Beaumont (2009) and Filippi (2013)*/
    for (k=0;k<dim*dim;k++)
    {
        L[k] *= SCALE;
    }
    /* really it is the Cholesky factorisation we need*/
    choldc(dim,L);
    return 0;
}

/** 
 * @brief proposal kernel sampler 
 * @details given kernel parameters L, and a particle \theta a new proposal 
 *          \theta* is generated.
 *
 * @param dim           dimensionality of parameter space
 * @param theta         particle from previous step
 * @param theta_prop    proposel particle update
 * @param L             proposal kernel parameters as generated by adapt_kernel
 */
int 
kernel_adapt(unsigned int dim,double *theta, double *theta_prop, double *L)
{
    double Z[dim];
    unsigned int tf;
    tf = 1;
    while (tf)
    {
        durngmvns(dim,theta,L,Z,theta_prop);
        tf = 0;
        for (int i=0;i<dim;i++)
        {
            tf = tf || (theta_prop[i] > abc_p.support[i+dim]) || (theta_prop[i] < abc_p.support[i]);
        }
    }
    return 0;
}

/** 
 * @brief proposal kernel probability density function
 * @details given kernel parameters L, a particle \theta and a new proposal 
 *          \theta* the density is evaluated pointwise.
 *
 * @param dim           dimensionality of parameter space
 * @param theta         particle from previous step
 * @param theta_prop    proposel particle update
 * @param L             proposal kernel parameters as generated by adapt_kernel
 *
 * @returns the density at point \theta* given \theta and L.
 */
double 
kernelPDF_adapt(unsigned int dim,double *theta, double *theta_prop, double *L)
{
    double x[dim];
    double y[dim];
    double z;
    double denom;
    unsigned int i;
    /* denom = sqrt(2pi)^d det(Sigma)*/ 
    denom = pow(2.0*M_PI,((double)dim));
    /* det(Sigma) = det(L)det(L^T) = det(L)^2 = prod_i^d L(i,i) */
    for (i=0;i<dim;i++)
    {
        denom *= L[i*dim + i]*L[i*dim+i];
    }

    for (i=0;i<dim;i++)
    {
        x[i] = theta_prop[i] - theta[i];
    }
    /*solve Ly = x */
    cholfs(dim, L, y, x);
    /* z = y^T y*/
    z = 0.0;
    for(i=0;i<dim;i++)
    {
        z += y[i]*y[i];
    }
    return exp(-0.5*z)/sqrt(denom);
}

/** 
 * @brief fixed proposal kernel sampler 
 * @details given a particle \theta a new proposal \theta* is generated.
 *
 * @param dim           dimensionality of parameter space
 * @param theta         particle from previous step
 * @param theta_prop    proposel particle update
 */
int 
kernel(unsigned int dim,double *theta, double *theta_prop)
{
    unsigned int i;
    for (i=0;i<dim;i++)
    {
        /* @todo must consider an adaptive scheme see Beaumont 2009 and Filippi 2013*/
        double prior_var = (abc_p.support[i+dim] - abc_p.support[i])*(abc_p.support[i+dim] - abc_p.support[i])/12.0;
        double sigma = sqrt(0.1*prior_var);
        theta_prop[i] = theta[i] + durngns(0.0,sigma);
        while (theta_prop[i] > abc_p.support[i+dim] || theta_prop[i] < abc_p.support[i])
        {
            theta_prop[i] = theta[i] + durngns(0.0,sigma);
        }
    }
    return 0;
}

/** 
 * @brief fixed proposal kernel probability density function
 * @details given a particle \theta and a new proposal \theta* the density is 
 *          evaluated pointwise.
 *
 * @param dim           dimensionality of parameter space
 * @param theta         particle from previous step
 * @param theta_prop    proposel particle update
 *
 * @returns the density at point \theta* given \theta.
 */
double 
kernelPDF(unsigned int dim,double *theta, double *theta_prop)
{
    unsigned int i;
    double p;
    p = 1;
    for (i=0;i<dim;i++)
    {
        
        double prior_var = (abc_p.support[i+dim] - abc_p.support[i])*(abc_p.support[i+dim] - abc_p.support[i])/12.0;
        double sigma2 = 0.1*prior_var;
        p *= exp(-((theta_prop[i] - theta[i])*(theta_prop[i] - theta[i]))/(2.0*sigma2))/sqrt(2.0*M_PI*sigma2);
    }
    return p; 
}


/**
 * @brief forwards simulation function of the data generation process.
 * @details i.i.d. samples form the Ornstein-Uhlenbeck SDE stationary distribution
 *
 * @param sim unused here for compatilibity with other codes
 * @param theta parameter values to use for simulation
 * @param data_s pointer to data structure to store simulated data
 */
int 
simulate_stat(void *sim,double * theta, Dataset * data_s)
{
    unsigned int N = 1000;
    double *X, *m_s;
    double D = theta[0];

    /*intial conditions*/
    m_s = (double*)data_s->fields[0].data_array;
    X = (double*)malloc(N*sizeof(double));
    /* sample the stationary distribution*/
    for (int i=0;i<N;i++) {
        X[i] = durngns(mu,sqrt(D/gam));
    }
    /*use Monte Carlo to estimate moments*/
    dmcint(N,1,X,NULL,NULL,m_s,m_s+1);
    approx_sim_counter++;
    free(X);
    return 0;
}

/**
 * @brief OU drift function
 */
void
OU_drift(double *X, unsigned int n, double *theta, unsigned int m, 
     double t, double *a) 
{
    for (int i=0;i<n;i++){
        a[i] = gam*(mu - X[i]);
    }
    return;
}

/**
 * @brief OU diffusion function
 */
void
OU_diffusion(double *X, unsigned int n, double *theta, unsigned int m, 
     double t, double *b) 
{
    double D = theta[0];
    for (int i=0;i<n;i++) {
        b[i] = sqrt(2.0*D);
    }
    return;
}


/**
 * @brief forwards simulation function of the data generation process.
 * @details Stochastic simulation of Ornstein-Uhlenbeck SDE
 *
 * @param sim unused here for compatilibity with other codes
 * @param theta parameter values to use for simulation
 * @param data_s pointer to data structure to store simulated data
 */
int 
simulate_SDE(void *sim,double * theta, Dataset * data_s)
{
    unsigned int N;
    N = 1000;
    double *X0;
    double *m_s;
    double *X;
    int *dims;
    dims = (int*)malloc(N*sizeof(int));
    X0 = (double*)malloc(N*sizeof(double));
    X = (double*)malloc(N*sizeof(double));
    
    /*intial conditions*/
    m_s = (double*)data_s->fields[0].data_array;

    /*initialise SDE at x0*/
    for (int i =0;i<N;i++)
    {
        X0[i] = x0;
        dims[i] = i; /*observe all sites*/
    }

    /*simulate N realisations of OU SDE using Euler-Maruyama*/
    daems(1,N,1,&T,theta,X0,&OU_drift,&OU_diffusion,N,dims,dt,X);
    /*use Monte Carlo to estimate moments*/
    dmcint(N,1,X,NULL,NULL,m_s,m_s+1);

    sim_counter++;
    free(X0);
    free(dims);
    free(X);
    return 0;
}

/**
 * @brief Import cell count data *.csv files
 *
 * @param filename name of *.csv file
 * @parm data pointer to data structure that will store the cell density profiles
 */
void
GenerateData(Dataset *data,double D)
{
    unsigned int N;
    N = 1000;
    double *X0;
    double *X;
    int *dims;
    dims = (int*)malloc(N*sizeof(int));
    X0 = (double*)malloc(N*sizeof(double));
    X = (double*)malloc(N*sizeof(double));
    double *m_d;

    /*allocate dataset*/
    data->numFields = 1; /* sum*/
    data->fields = (field *)malloc(data->numFields*sizeof(field));
    /*stores mean and var*/
    data->fields[0].numCols = 2;
    data->fields[0].numRows = 1;
    data->fields[0].type = REAL64_DATA;
    data->fields[0].numBytes = (data->fields[0].numCols)*sizeof(double);
    data->fields[0].data_array = malloc(data->fields[0].numBytes);
    memset(data->fields[0].data_array,0,data->fields[0].numBytes);
    
    m_d = (double*)data->fields[0].data_array;

/*initialise SDE at x0*/
    for (int i =0;i<N;i++)
    {
        X0[i] = x0;
        dims[i] = i; /*observe all sites*/
    }

   fprintf(stderr,"N=%d gam=%lg, mu=%lg, T=%lg D=%lg x0=%lg dt=%lg\n",N,gam,mu,T,D,x0,dt); 
    /*generate some data using Euler-Maruyama*/
    daems(1,N,1,&T,&D,X0,&OU_drift,&OU_diffusion,N,dims,dt,X);
    /*use Monte Carlo to estimate moments*/
    dmcint(N,1,X,NULL,NULL,m_d,m_d+1);


    return;
}


/**
 * Program entry point
 */
int 
main(int argc, char ** argv)
{
    double *theta,*weights;
    clock_t start_t, end_t;
    double time;
    int m = 1;
    double s;
    unsigned int N_particles;
    if (argc  < (8 + 2*m))
    {
        fprintf(stderr,"Usage: %s N eps T_smc T dt mu gamma x0 [sl1,sl2,...slm,su1,su2,sum] \n",argv[0]);
        exit(1);
    }
    else
    {
        /*set up ABC PC SMC params*/

        /* we need aabc_p = abc_p except for the simulation function*/
        /* number of particles*/
        N_particles = (unsigned int)atoi(argv[1]);
        
        abc_p.nacc = N_particles;
        aabc_p.nacc = N_particles;
        abc_p.nmax = 0;
        abc_p.nmax = abc_p.nmax;
        
        /* final target epsilon */
        abc_p.eps = (double)atof(argv[2]);
        aabc_p.eps = abc_p.eps ;
        
        /*set number of target distributions*/
        smc_p.T = (unsigned int)atoi(argv[3]);
        
        /* @todo sequence is geometric, but should also consider adaptive*/
        smc_p.eps_t = (double*)malloc((smc_p.T)*sizeof(double));
        smc_p.eps_t[smc_p.T-1] = abc_p.eps;
        for (int t=1;t<smc_p.T;t++)
        {
            smc_p.eps_t[smc_p.T-1-t] = 2.0*(smc_p.eps_t[smc_p.T - t]);
        }
        for (int t=0;t<smc_p.T;t++)
        {
            fprintf(stderr,"eps_t[%d] = %lf\n",t,smc_p.eps_t[t]);
        }
        
        /*parameters needed for the disrete model*/
        T = (double)atof(argv[4]);
        dt = (double)atof(argv[5]);
        mu = (double)atof(argv[6]);
        gam = (double)atof(argv[7]);
        x0 = (double)atof(argv[8]);
        
        /*fixed number of parameters*/
        abc_p.k = m; /* lambda, A and K*/
        aabc_p.k = abc_p.k;
        abc_p.support = (double *)malloc(2*abc_p.k*sizeof(double));
        aabc_p.support = (double *)malloc(2*aabc_p.k*sizeof(double));
        for (int k=0;k<abc_p.k*2;k++)
        {
            abc_p.support[k] = (double)atof(argv[9+k]);
        }
        memcpy(aabc_p.support, abc_p.support,2*abc_p.k*sizeof(double));
        abc_p.sim = NULL;
        aabc_p.sim = abc_p.sim;
        abc_p.rho = &rho;
        aabc_p.rho = abc_p.rho;
        abc_p.p = &prior;
        aabc_p. p = abc_p.p;
        abc_p.pd = &priorPDF;
        aabc_p.pd = abc_p.pd;
        abc_p.s = &simulate_SDE;
        aabc_p.s = &simulate_stat; 
        smc_p.E = (double)(N_particles);
        /*using adaptive kernels*/
        smc_p.q = &kernel;
        smc_p.qd = &kernelPDF;
        smc_p.q_params = (double *)malloc(abc_p.k*abc_p.k*sizeof(double));
        smc_p.q_adpt = &kernel_adapt;
        smc_p.qd_adpt = &kernelPDF_adapt;
        smc_p.adpt = &adapt_kernel;
    }
    /*initialise our RNG library*/
    SSAL_Initialise(argc,argv);
    /* import data*/
    GenerateData(&data,10.0);

    /* Sanity check: print data summary */
    {
        double *m_d;
        int ns;
        ns = data.fields[0].numCols;
        m_d = (double*)data.fields[0].data_array;
        for (int j=0;j<ns;j++)
        {
            fprintf(stderr,"S_%d = %lg\n",j,m_d[j]);
        }
   //     int nt;
   //     nt = data.fields[0].numCols;
   //     T = (double *)data.fields[0].data_array;
   //     c0 = *((double *)data.fields[1].data_array);
   //     C_dat = (double *)data.fields[2].data_array;
   //     fprintf(stderr, "nt = %d, C0 = %lg\n",nt, c0);
   //     fprintf(stderr,"T C\n");
   //     for (int j=0;j<nt;j++)
   //     {
   //         fprintf(stderr,"%lg %lg\n",T[j],C_dat[j]);
   //     }
   //     
   //     /* check that any sample of with be K > c0*/
   //     if (abc_p.support[2] <= c0)
   //     {
   //         fprintf(stderr,"Invalid Support for K, ensure C(0) < K for any K ~ p(K).\n C(0) = %lg, p(K) = U(%lg,%lg).\n",c0,abc_p.support[2],abc_p.support[5]);
   //         exit(1);
   //     }
   }
   // 

    /*allocate output array*/
    theta = (double *)malloc(N_particles*abc_p.k*sizeof(double));
    weights = (double *)malloc(N_particles*sizeof(double));

    approx_sim_counter = 0; /*for performance metric*/
    sim_counter = 0; /*for performance metric*/

    /*generate posterior samples using ABC-SMC with Preconditiong*/
    start_t = clock();
    dabcapcsmc(aabc_p,abc_p,smc_p,0,&data,theta,weights,NULL);
    end_t = clock();

    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    
    /*write output*/
    fprintf(stdout,"\"Particle\",\"D\",\"weights\",\"SEC\",\"NACC\",\"NASIMS\",\"NSIMS\"\n");
    for (unsigned int i=0;i<N_particles;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%u,%u,%u\n",i,theta[i],weights[i],time,abc_p.nacc,approx_sim_counter,sim_counter);
    }

    return 0;
}

