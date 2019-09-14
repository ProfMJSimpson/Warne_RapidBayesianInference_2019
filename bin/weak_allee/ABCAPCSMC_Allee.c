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

/* ABC Preconditioned SMC inference using the generalised logistic growth model 
 * for the preconditioner and final target of the the hexagonal lattice model 
 * with a generalised proliferation function of Jin et al. (2016) Phys. Biol. 
 * using cell density profiles.
 */
ABC_Parameters aabc_p; /* for the continuum model*/
ABC_Parameters abc_p; /* for the discrete model*/
SMC_Parameters smc_p;


/* data set structure */
Dataset data;

/* fixed parameters for the discrete model*/
double Pm;
int I;
int J;

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
    double *C,*C_s;
    double d;
    unsigned int nt;
    
    C_s = (double *)data_s->fields[2].data_array;
    C = (double *)data->fields[2].data_array;

    nt = (unsigned int)data->fields[0].numCols;

    d = 0;
    for (unsigned int i=0;i<nt;i++)
    {
        d += (C[i] - C_s[i])*(C[i] - C_s[i]);
    }
    return sqrt(d);
}

/**
 * @brief prior sampler 
 * @details Generates a number of i.i.d. samples from p(\lambda,A,K)
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
        theta[i*dim] = durngus(support[0],support[3]); // lambda 
        theta[i*dim + 2] = durngus(support[2],support[5]); // K
        theta[i*dim + 1] = durngus(support[1],theta[i*dim+2]); // A <= K
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
    tf = tf && (theta[1] <= theta[2]); // ensure A <= K
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
 * @brief Implements S(c) = \lambda*c*f(c) for the weak Allee effect continuum
 *        model.
 * @param c         pointer to a population density value
 * @param n         number of density points
 * @param params    model parameters
 * @param f_r       pointer to store S(c)
 */
void 
wae(double *c, unsigned int n, double *params,unsigned int m, double t,double *f_r)
{
    double lambda, A, K;
    lambda = params[0]; 
    A = params[1];
    K = params[2];
    f_r[0] = lambda*c[0]*(1.0 - c[0]/K)*((c[0] + A)/K);
}

/**
 * @brief neighbourhood function for hexagonal lattice with periodic boundaries
 *
 * @param d         dimensionality of lattice
 * @param N         array of number of sites along each lattice axis
 * @param ind       index of a lattice site
 * @param Nhsize    pointer to store the number of sites that are neighbours to 
 *                  site ind
 * @param hex       array of size indices that are neighbours to site ind
 */
void
hexNh(int d, int *N, int ind, int *Nhsize, int *hex)
{
    int Nx = N[0];
    int Ny = N[1];
    int offx_oddy[6] = {0,-1,0,1,1,1};
    int offx_eveny[6] = {-1,-1,-1,0,1,0};
    int offy[6] = {-1,0,1,1,0,-1};
    
    /* get 2d index*/
    int j = ind / Nx;
    int i = ind % Nx;

    *Nhsize = 6;
    /*implementing periodic boundary*/
    if (j%2 == 0)
    {
        for (int k=0;k<*Nhsize;k++)
        {   
            hex[k] = ((j+offy[k]+Ny)%Ny)*Nx + (i+offx_eveny[k]+Nx)%Nx;
        }
    }
    else
    {
        for (int k=0;k<*Nhsize;k++)
        {   
            hex[k] = ((j+offy[k]+Ny)%Ny)*Nx + (i+offx_oddy[k]+Nx)%Nx;
        }
    }
}

/** 
 * @brief motility function
 *
 * @param Cavg      local density
 * @param params     model parameters
 *
 * @returns the probability of a successful motility event
 */
double 
g(double Cavg, double *params)
{
    return Pm; /*leads to linear diffusion in continuum limit*/
}

/** 
 * @brief crowding function
 *
 * @param Cavg      local density
 * @param params     model parameters
 *
 * @returns the probability of a successful proliferation event
 */
double 
f(double Cavg, double *params)
{
    double Pp = params[0];
    double A = params[1];
    double K = params[2];
    /* leads to the Weak Allee effect in the continuum limit 
     * dC/dt = \lambda C(1-C/K)((C + A)/K)
     */
    return Pp*(1.0 - Cavg/K)*((Cavg + A)/K); 
}

/**
 * @brief forwards simulation function of the data generation process.
 * @details Numerically solves the ODE for the continuum model
 *
 * @param sim unused here for compatilibity with other codes
 * @param theta parameter values to use for simulation
 * @param data_s pointer to data structure to store simulated data
 */
int 
simulate_ODE(void *sim,double * theta, Dataset * data_s)
{
    unsigned int m,nt;
    double *C_r;
    double C0;;
    double *T; 
    unsigned int d = 0;
    double h;

    /*intial conditions*/
    T = (double*)data.fields[0].data_array;
    nt = data.fields[0].numCols;
    C0 = *((double*)data.fields[1].data_array);
    C_r = (double*)data_s->fields[2].data_array;

    h = T[nt-1]/10000.0;
    m = 3; /*params are lambda, beta, and gamma */
    
    /*solve ODE with RKF45 method*/
    drkf45s(m,1,nt,T,theta,&C0,&wae,1,&d,h,1e-8,C_r);

    approx_sim_counter++;
    return 0;
}

/**
 * @brief forwards simulation function of the data generation process.
 * @details Stochastic simulation of lattice-based discrete random walk model
 *
 * @param sim unused here for compatilibity with other codes
 * @param theta parameter values to use for simulation
 * @param data_s pointer to data structure to store simulated data
 */
int 
simulate_DRW(void *sim,double * theta, Dataset * data_s)
{
    unsigned int m,nt;
    double *C_r,*C0, *C;
    double c0;;
    double *T; 
    int *dims;
    int N[2];

    N[0] = I;
    N[1] = J;

    /*intial conditions*/
    T = (double*)data.fields[0].data_array;
    nt = data.fields[0].numCols;
    c0 = *((double*)data.fields[1].data_array);
    C = (double*)data_s->fields[2].data_array;

    m = 3; /*params are lambda, A, and K */
    
    /*initialise lattice and initial occupacies*/
    C0 = (SSAL_real_t*) malloc(N[0]*N[1]*sizeof(double));
    dims = (int*)malloc(N[0]*N[1]*sizeof(int));
    for (int i =0;i<N[0]*N[1];i++)
    {
        SSAL_real_t u;
        u = durngus(0,1);
        C0[i] =  (SSAL_real_t)(u <= c0);
        dims[i] = i; /*observe all sites*/
    }

    C_r = (SSAL_real_t*)malloc(N[0]*N[1]*nt*sizeof(double));
    memset((void*)C_r,0,N[0]*N[1]*nt*sizeof(double));
    memset((void*)C,0,nt*sizeof(double));
    /*simulate stochastic realisations*/
    dalrws(2,N,&hexNh,6,nt,T,theta,C0,&g,&f,N[0]*N[1],dims,1.0,C_r);
    
    /* take spatial average*/
    for (int i=0;i<N[0]*N[1];i++)
    {
        for (int ti=0;ti<nt;ti++)
        {
            C[ti] += C_r[i*nt + ti];        
        }
    }
    for (int ti=0;ti<nt;ti++)
    {
        C[ti] /= (double)(N[0]*N[1]);        
    }
    sim_counter++;
    free(C0);
    free(dims);
    free(C_r);
    return 0;
}

/**
 * @brief Import cell count data *.csv files
 *
 * @param filename name of *.csv file
 * @parm data pointer to data structure that will store the cell density profiles
 */
void
ImportCellDensityData(char * filename, Dataset *data)
{
    FILE *fp;
    int numRows, numCols;
    double *buffer;
    char temp[125];
    double *T, *C, *C0;
    /*open file*/
    fp = fopen(filename,"r");
    /* count rows and columns*/
    numCols = 0;
    numRows = 0;
    while(!feof(fp))
    {
        char chr = fgetc(fp);
        numCols += (chr == ' ' && numRows == 0);
        numRows += (chr == '\n');
    }
    /*back to start of file fo parsing*/
    rewind(fp);
    /*buffer for reading*/
    buffer = (double *)malloc(numRows*numCols*sizeof(double));
    for (int i=0;i<numRows;i++)
    {
        fscanf(fp,"%s",temp);
        for (int j=0;j<numCols;j++)
        {
            fscanf(fp,"%lf",buffer + (i*numCols + j));
        }
        fscanf(fp,"\n");
    }
    /*close the file now*/
    fclose(fp);

    /*allocate dataset*/
    data->numFields = 3; /* time-axis, initial density, density*/
    data->fields = (field *)malloc(data->numFields*sizeof(field));
    /*time axis*/
    data->fields[0].numCols = numCols-1;
    data->fields[0].numRows = 1;
    data->fields[0].type = REAL64_DATA;
    data->fields[0].numBytes = (data->fields[0].numCols)*sizeof(double);
    /* initial density */
    data->fields[1].numCols = 1;
    data->fields[1].numRows = 1;
    data->fields[1].type = REAL64_DATA;
    data->fields[1].numBytes = (data->fields[1].numCols)*sizeof(double);
    /*density profile */
    data->fields[2].numCols = numCols-1;
    data->fields[2].numRows = 1;
    data->fields[2].type = REAL64_DATA;
    data->fields[2].numBytes = (data->fields[2].numCols)*sizeof(double);

    for (int i=0;i<data->numFields;i++)
    {
        data->fields[i].data_array = malloc(data->fields[i].numBytes);
        memset(data->fields[i].data_array,0,data->fields[i].numBytes);
    }

    /* copy data to data structure*/
    T = (double *)data->fields[0].data_array;
    C0 = ((double*)data->fields[1].data_array);
    C = (double*)data->fields[2].data_array;
    
    for (int j=1;j<numCols;j++)
    {
        T[j-1] = buffer[j];
    }

    for (int i=1;i<numRows;i++)
    {
        C0[0] += buffer[i*numCols];
        for (int j=1;j<numCols;j++)
        {
            C[j-1] += buffer[i*numCols + j];
        }
    }
    C0[0] /= (double)(numRows-1);
    for (int j=0;j<numCols-1;j++)
    {
        C[j] /= (double)(numRows-1);
    }
    free(buffer);
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
    char * filename;
    int m = 3;
    double s;
    unsigned int N_particles;
    if (argc  < (8 + 2*m))
    {
        fprintf(stderr,"Usage: %s filename N eps T Pm I J [sl1,sl2,...slm,su1,su2,sum] \n",argv[0]);
        exit(1);
    }
    else
    {
        filename = argv[1];
        /*set up ABC PC SMC params*/

        /* we need aabc_p = abc_p except for the simulation function*/
        /* number of particles*/
        N_particles = (unsigned int)atoi(argv[2]);
        
        abc_p.nacc = N_particles;
        aabc_p.nacc = N_particles;
        abc_p.nmax = 0;
        abc_p.nmax = abc_p.nmax;
        
        /* final target epsilon */
        abc_p.eps = (double)atof(argv[3]);
        aabc_p.eps = abc_p.eps ;
        
        /*set number of target distributions*/
        smc_p.T = (unsigned int)atoi(argv[4]);
        
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
        Pm = (double)atof(argv[5]);
        I = (int)atoi(argv[6]);
        J = (int)atoi(argv[7]);
        
        /*fixed number of parameters*/
        abc_p.k = m; /* lambda, A and K*/
        aabc_p.k = abc_p.k;
        abc_p.support = (double *)malloc(2*abc_p.k*sizeof(double));
        aabc_p.support = (double *)malloc(2*aabc_p.k*sizeof(double));
        for (int k=0;k<abc_p.k*2;k++)
        {
            abc_p.support[k] = (double)atof(argv[8+k]);
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
        abc_p.s = &simulate_DRW;
        aabc_p.s = &simulate_ODE; 
        smc_p.E = (double)(N_particles);
        /*using adaptive kernels*/
        smc_p.q = &kernel;
        smc_p.qd = &kernelPDF;
        smc_p.q_params = (double *)malloc(abc_p.k*abc_p.k*sizeof(double));
        smc_p.q_adpt = &kernel_adapt;
        smc_p.qd_adpt = &kernelPDF_adapt;
        smc_p.adpt = &adapt_kernel;
    }
    
    /* import data*/
    ImportCellDensityData(filename,&data);

    /* Sanity check: print data summary */
    {
        double *T, *C_dat, c0;
        int nt;
        nt = data.fields[0].numCols;
        T = (double *)data.fields[0].data_array;
        c0 = *((double *)data.fields[1].data_array);
        C_dat = (double *)data.fields[2].data_array;
        fprintf(stderr, "nt = %d, C0 = %lg\n",nt, c0);
        fprintf(stderr,"T C\n");
        for (int j=0;j<nt;j++)
        {
            fprintf(stderr,"%lg %lg\n",T[j],C_dat[j]);
        }
        
        /* check that any sample of with be K > c0*/
        if (abc_p.support[2] <= c0)
        {
            fprintf(stderr,"Invalid Support for K, ensure C(0) < K for any K ~ p(K).\n C(0) = %lg, p(K) = U(%lg,%lg).\n",c0,abc_p.support[2],abc_p.support[5]);
            exit(1);
        }
    }
    
    /*initialise our RNG library*/
    SSAL_Initialise(argc,argv);

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
    fprintf(stdout,"\"Particle\",\"lambda\",\"A\",\"K\",\"weights\",\"SEC\",\"NACC\",\"NASIMS\",\"NSIMS\"\n");
    for (unsigned int i=0;i<N_particles;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%g,%g,%u,%u,%u\n",i,theta[i*abc_p.k],theta[i*abc_p.k+1],theta[i*abc_p.k + 2],weights[i],time,abc_p.nacc,approx_sim_counter,sim_counter);
    }

    return 0;
}

