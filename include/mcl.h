/* MCL: Monte Carlo Library 
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

#ifndef __LIBABC_H_
#define __LIBABC_H_

#include "SSAL.h"
#include <math.h>
#include <stdint.h>
#include <time.h>

#define MCL_MAX_NAME_SIZE 128
enum DataType_enum
{
    REAL32_DATA,
    REAL64_DATA,
    INT32_DATA,
    INT64_DATA,
    NAT32_DATA,
    NAT64_DATA,
    TEXT_DATA,
    PTR_DATA
};
typedef enum DataType_enum DataType;

#if defined(__CHECKPOINT__)
#define CHECKPOINT_FILENAME "sim.chkpnt"
#endif
extern unsigned int sim_counter; 
extern unsigned int approx_sim_counter; 

struct field_struct 
{
    /**field name*/
    char name[MCL_MAX_NAME_SIZE];    
    size_t numRows; 
    size_t numCols;
    size_t numBytes;
    DataType type;
    void *data_array;

};
typedef struct field_struct field;

/* generic dataset structure*/
struct Dataset_struct 
{
    size_t numFields;
    field *fields;
};
typedef struct Dataset_struct Dataset;

struct ABC_Parameters_struct 
{ 
    /**total max samples before termination*/
    unsigned int nmax; 
    /**If nmax is used, then use percentile to INFER eps*/
    SSAL_real_t percentile;
    /**acceptance threshold*/   
    SSAL_real_t eps; 
    /**number of acceptances required*/ 
    unsigned int nacc; 
    /**dimensionality of parameter space*/  
    unsigned int k; 
    /**prior sampling region*/ 
    SSAL_real_t *support; 
    /**Model simulation object */ 
    void *sim; 
    /**distance function*/ 
    SSAL_real_t (*rho)(Dataset *,Dataset *); 
    /**prior sampler function*/ 
    int (*p)(unsigned int,unsigned int,SSAL_real_t *,SSAL_real_t *); 
    /** Prior Distribution function; not required for all methods*/ 
    SSAL_real_t (*pd)(unsigned int, SSAL_real_t*); 
    /**model simulation function*/ 
    int (*s)(void *, SSAL_real_t* ,Dataset *); 
};
typedef struct ABC_Parameters_struct ABC_Parameters;

struct SMC_Parameters_struct 
{
    /**number of steps*/
    unsigned int T;
    /**ABC thresholds*/
    SSAL_real_t *eps_t; 
    /**Effective sample size (ESS) threshold for re-sampling.*/
    SSAL_real_t E; 
    /**transition kernel sampler*/ 
    int (*q)(unsigned int, SSAL_real_t*, SSAL_real_t*); 
    /**transition kernel Distribution function*/
    SSAL_real_t (*qd)(unsigned int, SSAL_real_t*, SSAL_real_t*); 
    /**kernel parameters for adaptive proposals*/
    SSAL_real_t *q_params;
    /**adaptive transition kernel sampler*/ 
    int (*q_adpt)(unsigned int, SSAL_real_t*, SSAL_real_t*, SSAL_real_t*); 
    /**adaptive transition kernel Distribution function*/
    SSAL_real_t (*qd_adpt)(unsigned int, SSAL_real_t*, SSAL_real_t*, SSAL_real_t*);
    /**proposal adaptation method*/
    int (*adpt)(unsigned int, unsigned int, SSAL_real_t*, SSAL_real_t*, SSAL_real_t*);
};
typedef struct SMC_Parameters_struct SMC_Parameters;

/*ABC samplers*/
int 
dabcrs(ABC_Parameters, Dataset *, double *,double *);

int 
dabcasmc(ABC_Parameters, SMC_Parameters, Dataset *, double *, double *,double *);

int
dabcapcsmc(ABC_Parameters, ABC_Parameters, SMC_Parameters, int, Dataset *,
          double *, double*, double *);
int 
dabcammsmc(ABC_Parameters, ABC_Parameters, SMC_Parameters, Dataset *, double *, 
           double *, double *, double *, double *, double *);

/*Monte Carlo integration functions*/
int 
dmcint(unsigned int, unsigned int, double *, double*, 
       double (*)(int , double *,double *), double *, double *); 

int 
dmcintd(unsigned int, unsigned int, double *, double *, double *,
        double (*)(int , double *, double*), double *, double *); 

int 
dmcintv(unsigned int, unsigned int, double *,double *, double *);

int 
dmcintcv(unsigned int, unsigned int, double *,double *, double *);

/*data manipulation*/
Dataset *
copyDataset(Dataset *);

/*some utilies*/
int
choldc(unsigned int, double *);

int
cholfs(unsigned int, double *, double *, double *);

int
mvprod(unsigned int, double *, double *, double *);
#endif
