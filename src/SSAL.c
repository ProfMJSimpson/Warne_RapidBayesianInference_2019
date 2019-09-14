/* SSAL: Stochastic Simulation Algorithm Library
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
/**
 * @file SSAL.c
 * @brief A fast implementation of standard Stochastic Simulation algorithms.
 * @details This library provides optimised sequential and parallel exact and 
 * approximate stochastic simulation algorithms. The implementations provide an
 * easy to use application program interface (API) using standard C structures.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Science and Engineering Faculty
 * @author Queensland University of Technology
 *
 */

#include "SSAL.h"

/**
 * @brief Handle error reporting
 * @param errCode The Error code 
 * @param funcName The Function namee in which the error occurred
 * @param lineNum The line number at which the error occurred
 * @param fatal Flag indicating the error was fatal and the progam should be 
 *              aborted
 * @param warning Indicates the error is to be interpreted as a warning
 * @param msg A custom message which can be appended to the default
 */
void 
SSAL_HandleError(int errCode, char * funcName, int lineNum, unsigned char fatal, 
                 unsigned char warning, char *msg)
{
    char * errName;
    char * defaultMsg;
    char * errType;
    switch (errCode)
    {
        case SSAL_UNKNOWN_OPTION_ERROR:
            errName = "SSAL_UNKNOWN_OPTION_ERROR";
            defaultMsg = "A configuration option was not recognised.";
            break;
        default:
            errName = "SSAL_UNKNOWN_ERROR";
            defaultMsg = "No idea what happened then...";
            break;
    }

    if (warning)
    {
        errType = "ERROR";   
    }
    else
    {
        errType = "WARNING";
    }
        
    if (msg != NULL)
    {
        fprintf(stderr,"[%s %s]: %s\n %s\n Error Code [%d] Line [%d]\n",
            errType,errName,defaultMsg,msg,errCode,lineNum);    
    }
    else
    {
        fprintf(stderr,"[%s %s]: %s\n Error Code [%d] Line [%d]\n",
            errType,errName,defaultMsg,errCode,lineNum);    
    }
    

    if (fatal)
    {
        fprintf(stderr,"Fatal Error Aborting!\n");
        exit(1);
    }
    return;
}



/**
 * @brief Initilialises the SSAL library
 * @detail Initialised the computational backend based on compilation
 * options applied. Some user customisations are also available.
 *
 * @param argc the number of input args
 * @param argv an array of args
 * @todo Will add more options here as needed.
 * @retVal SSAL_SUCCESS if initialisation was successful
 */
int 
SSAL_Initialise(int argc,char **argv)
{
    int i;
    unsigned char infoflag;
    unsigned int seed;
    seed = 1337;
    infoflag = 0;
    /*search for config options*/
    for (i=1;i<argc;i++)
    {
        if (!strcmp("--SSALOPT",argv[i]))
        {
            i++;
            break;
        }
    }
    

    /*set defaults*/
    /*parse args*/
    for (;i<argc;i++)
    {
        if (!strcmp("--seed",argv[i]))
        {
            seed = (unsigned int)atoi(argv[++i]);          
        }
        else if (!strcmp("--info",argv[i]))
        {
            infoflag = 1;
        }
        else
        {
            SSAL_HandleError(SSAL_UNKNOWN_OPTION_ERROR,"SSAL_Initialise",
                             __LINE__,0,1,argv[i]);
        }
    }

    /*initialise RNGS */
    #ifdef __SERIAL__
       suarngs(seed,NULL,NULL);
    #endif
    if (infoflag)
    {
        fprintf(stderr,"RNG Seed: %d\n",seed);
        fprintf(stderr,"RNG RAND_MAX: %d\n",RAND_MAX);
    }
    return SSAL_SUCCESS;
}

