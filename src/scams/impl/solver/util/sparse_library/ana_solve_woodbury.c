/*****************************************************************************

    Copyright 2010 - 2013
    Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V.


   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

 *****************************************************************************/

/*****************************************************************************

 ana_solve_woodbury.c - description

 Original Author: Christiane Reuther Fraunhofer IIS/EAS Dresden

 Created on: 20.09.2011

 SVN Version       :  $Revision: 1754 $
 SVN last checkin  :  $Date: 2014-05-26 11:18:44 +0200 (Mon, 26 May 2014) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: ana_solve_woodbury.c 1754 2014-05-26 09:18:44Z karsten $

 *****************************************************************************/

/**
 * @file 	ana_solve_woodbury.c
 * @author	Christiane Reuther
 * @date	November 07, 2012
 * @brief	Source-file to define method <i>ana_solve_woodbury</i>
 */

/*****************************************************************************/


#include "ana_solv_data.h"
#include "ma_typedef.h"
#include "ma_util.h"
#include "ma_sparse.h"
#include "ma_woodbury.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 * The method <i>ana_solve_woodbury</i> computes the solution of a linear system
 * of equations by calling the function <i>MA_LowRankModifications</i>:
 * <ul>
 * <li> In case of <i>sdata->cur_algorithm = EULER</i>: \f$x\f$ is the solution to
 * \f$\left(Z_{euler} + V\,W\right) x = W_{euler}\, x_{last} - q\f$
 * <li> In case of <i>sdata->cur_algorithm = TRAPEZ</i>: \f$x\f$ is the solution to
 * \f$\left(Z_{trapez} + V\,W\right) x = W_{trapez}\, x_{last} + A\,xp - q\f$
 * </ul>
 *
 * In case of variable step size (<i>sdata->variable_step_size</i> = 1): replace
 * the sparse matrices \f$W_{euler}\f$ and \f$W_{trapez}\f$ by
 * <ul>
 * <li>	\f$W_{euler} = \frac1{h_{temp}} A\f$,
 * <li>	\f$W_{trapez} = \frac2{h_{temp}} A\f$,
 * </ul>
 * where \f$h_{temp} = \f$ <i>sdata->h_temp</i> is the current step size.
 *
 *  Implementation details of the function <i>ana_solv</i>: <br>
 *	<ul>
 *	<li> The vector \f$ x \f$ deals as input as solution vector \f$x_{last}\f$
 *	of the previous time step.
 *	<li> The vector \f$ x \f$ deals as output as solution vector of the current
 *	time step.
 *	<li> The vector \f$xp\f$ is computed by \f$xp = \frac1{h}(x - x_{last})\f$
 *	after solving the linear system of equation, \f$xp\f$ is used as input at
 *	the next time step.
 *	<li>  The parameters <i>k, mode, n, s</i> are used within the function
 *  <i>MA_LowRankModifications</i> to define the matrices \f$V\f$ and \f$W\f$.
 *	<li> In case of <i>first_call</i> = 1 allocate memory for
 *	<i>sdata->aux_woodbury</i> which is used in the function
 *	<i>MA_LowRankModifications</i>.
 *	</ul>
 *
 *	@return
 *  <ul><li>    0 - okay
 *  <li>        1 - reallocation of NULL pointer
 *  <li>        2 - not enough memory
 *  <li>		3 - dimension erroneous
 *  <li>		4 - matrix singular
 *  <li>		6 - no sparse matrix
 *  </ul>
 */

int ana_solve_woodbury(
          double* q,
          double* x,
          sca_solv_data* sdata,
          long k,
          int* mode,
          long* n,
          double* s,
          int first_call
      )
{
	double hinv;
	double *r1, *r2;
	unsigned long size;
	unsigned long i;
	double *xp, *x_last;

	struct sparse* W_temp = NULL;

	int err = 0;
	long dim;

	r1    = sdata->r1;
	r2    = sdata->r2;
	size  =sdata->size;
	xp    =sdata->xp;
	x_last=sdata->x_last;

   memcpy(x_last, x, size*sizeof(double));

   dim = 2*sdata->size*sdata->size;


   if (first_call != 0)
   {
  	   sdata->aux_woodbury =
  			   (double *) calloc((unsigned) dim, (unsigned) sizeof(double));

  	   if (sdata->aux_woodbury == NULL)
  		   return 2;
   }
   else
   {
	   memset(sdata->aux_woodbury,0,dim*sizeof(double));
   }

   /*check for to small stepsize. Abort if divisor is zero. */
   if(sdata->variable_step_size == 0)   /* fixed stepsize */
   {
	   if(fabs(sdata->h) < 1e-300 )
	   {
		   return -1;
		  /*printf("%s \n", "Division by zero in ana_solve_woodbury.c stepsize(1/h) to small! Aborting...");
		  //exit(-1); */
	   }
   }else	/* variable stepsize */
   {
	   if(fabs(sdata->h_temp) < 1e-300 )
	   {
		   return -1;
		   /*printf("%s \n", "Division by zero in ana_solve_woodbury.c variable stepsize(1/h_temp) to small! Aborting...");
		   //exit(-1); */
	   }
   }

   if(sdata->cur_algorithm==EULER)
   {  /******* Euler backward ********/
	  if (sdata->variable_step_size == 0)
	  {
		  hinv  = 1/sdata->h;

		  err = MA_ProductSparseVector(sdata->sW_euler, x, r1);
		  if (err)
			  return err;
	  }
	  else /* variable stepsize is set */
	  {
		  hinv  = 1/sdata->h_temp;

	      W_temp = (struct sparse*) calloc(1,(unsigned) sizeof (struct sparse));

	      MA_GenerateProductValueSparse(W_temp, sdata->A, hinv);
	   	  MA_ProductSparseVector(W_temp, x, r1);

	   	  free(W_temp);
	  }

      for(i=0;i<size;++i) r1[i] -= q[i];

      MA_LequSparseSolut(sdata->sZ_euler, sdata->code_euler, r1, x);
      sdata->cur_algorithm = EULER;

	  err = MA_LowRankModifications(sdata, x, k, mode, n, s);
   	  sdata->variable_step_size = 0;

      if(err)
      {
    	 sdata->critical_column=sdata->code_trapez->critical_column;
     	 sdata->critical_row=sdata->code_trapez->critical_line;
     	 return(err);
      }

      for(i=0;i<size;++i) xp[i]=(x[i]-x_last[i])*hinv; /*derivation*/

      sdata->cur_algorithm=TRAPEZ;
   }
   else
   {  /*** trapezoidal method ***/
	  if (sdata->variable_step_size == 0)
	  {
		  hinv=2.0/sdata->h;

		  MA_ProductSparseVector(sdata->sW_trapez, x, r1);     /* W*x(i-1) */
	  }
	  else /* variable stepsize is set */
	  {
		  hinv=2.0/sdata->h_temp;

	      W_temp = (struct sparse*) calloc(1,(unsigned) sizeof (struct sparse));

	      MA_GenerateProductValueSparse(W_temp, sdata->A, hinv);
	   	  MA_ProductSparseVector(W_temp, x, r1);

	   	  free(W_temp);
	  }

      MA_ProductSparseVector(sdata->A, xp, r2);            /* A*xp(i-1) */

      for(i=0;i<size;++i)
    	  r1[i] += r2[i] - q[i];     /* W*x(i-1) + A*xp(i-1) - q(i) */

      MA_LequSparseSolut(sdata->sZ_trapez, sdata->code_trapez, r1, x);

	  err = MA_LowRankModifications(sdata, x, k, mode, n, s);
   	  sdata->variable_step_size = 0;

      if(err)
      {
     	 sdata->critical_column=sdata->code_trapez->critical_column;
     	 sdata->critical_row=sdata->code_trapez->critical_line;
     	 return(err);
      }

      for(i=0;i<size;++i) xp[i]=hinv*(x[i]-x_last[i])-xp[i]; /*new derivation*/
   }


   return err;
}
