/*****************************************************************************

    Copyright 2010-2014
    Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V.

    Copyright 2015-2016
    COSEDA Technologies GmbH

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

 ana_utilities.c - description

 Original Author: Christiane Reuther Fraunhofer IIS/EAS Dresden

 Created on: 27.07.2012

 SVN Version       :  $Revision: 1927 $
 SVN last checkin  :  $Date: 2016-02-26 12:32:21 +0100 (Fri, 26 Feb 2016) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: ana_utilities.c 1927 2016-02-26 11:32:21Z karsten $

 *****************************************************************************/

/**
 * @file 	ana_utilities.c
 * @author	Christiane Reuther
 * @date	November 07, 2012
 * @brief	Source-file to define methods <i>ana_get_error_position</i>,
 * <i>ana_get_algorithm</i>, <i>ana_set_algorithm</i>, <i>ana_get_dimension</i>,
 * <i>ana_get_pivot_flop</i>, <i>ana_get_dec_flop</i>, <i>ana_get_sol_flop</i>,
 * <i>ana_set_variable_step_size</i>, <i>ana_store_check_point</i>,
 * <i>ana_restore_check_point</i>,  <i>ana_store_intermediate_check_point</i>
 * and <i>ana_restore_intermediate_check_point</i>
 */

/*****************************************************************************/


#include "ana_solv_data.h"
#include "ma_typedef.h"
#include "ma_util.h"
#include "ma_sparse.h"
#include <string.h>
#include <stdlib.h>


/****************************************/

void ana_get_error_position(sca_solv_data* data,long* row,long* column)
{
	(*row)=data->critical_row;
	(*column)=data->critical_column;
}

/****************************************/

int ana_get_algorithm(sca_solv_data* data)
{
	if (data->cur_algorithm == TRAPEZ) return 2;
	else return 1;
}

/****************************************/

void ana_set_algorithm(sca_solv_data* data, int alg)
{
	if (alg == 0) data->algorithm = EULER;
	else if (alg == 1) data->algorithm = TRAPEZ;
}

/****************************************/

long int ana_get_dimension(sca_solv_data* data)
{
	return data->size;
}

/****************************************/

long int ana_get_pivot_flop(sca_solv_data* data, int alg)
{
	if (alg == 1) return data->code_euler->pivot_flop;
	else return data->code_trapez->pivot_flop;
}

/****************************************/

long int ana_get_dec_flop(sca_solv_data* data, int alg)
{
	if (alg == 1) return data->code_euler->dec_flop;
	else return data->code_trapez->dec_flop;
}

/****************************************/

long int ana_get_sol_flop(sca_solv_data* data, int alg)
{
	if (alg == 1) return data->code_euler->sol_flop;
	else return data->code_trapez->sol_flop;
}

/****************************************/

void ana_set_variable_step_size(sca_solv_data* data, double h_new,
		double h_diff)
{
	data->variable_step_size = 1;
	data->h_temp = h_new;
	data->h_diff = h_diff;
}

/****************************************/


int ana_store_solver_check_point(
		  sca_solv_data* sdata,		            /**< internal solver data */
		  double* x,					        /**< solution vector */
		  sca_solv_checkpoint_data** cp_data    /**< internal check point data */
		  )
{
	unsigned long size;
	unsigned long dsize;

	if(sdata==NULL) return 1;

	size=sdata->size;
	dsize=size*sizeof(double);

	if(*cp_data==NULL)
	{
		(*cp_data)=(sca_solv_checkpoint_data*)malloc(sizeof(sca_solv_checkpoint_data));
		if((*cp_data)==NULL) return 2;
		(*cp_data)->size=size;

		if(size>0)
		{
			(*cp_data)->x= (double*)malloc(dsize);
			(*cp_data)->xp=(double*)malloc(dsize);
		}
		else
		{
			(*cp_data)->x=NULL;
			(*cp_data)->xp=NULL;
			(*cp_data)->size=0;
		}
	}
	else
	{
		if(size!=((*cp_data)->size))
		{
			if((*cp_data)->x!=NULL)
			{
				double* x=(double*)realloc((void*)(*cp_data)->x,dsize);
				if(x==NULL) free((*cp_data)->x);
				(*cp_data)->x=x;
			}
			else                    (*cp_data)->x=(double*)malloc(dsize);
			if((*cp_data)->xp!=NULL)
			{
				double* xp=(double*)realloc((void*)(*cp_data)->xp,dsize);
				if(xp==NULL) free((*cp_data)->xp);
				(*cp_data)->xp=xp;
			}
			else                     (*cp_data)->xp=(double*)malloc(dsize);

			(*cp_data)->size=size;
		}
	}

	if(size==0) return 0;

	/*can't allocate enough memory*/
	if(((*cp_data)->x==NULL) || ((*cp_data)->xp==NULL))
	{
		if((*cp_data)->x!=NULL) free((*cp_data)->x);
		(*cp_data)->x=NULL;

		if((*cp_data)->xp!=NULL) free((*cp_data)->xp);
		(*cp_data)->xp=NULL;

		(*cp_data)->size=0;

		return 2;
	}

	/*store x and xp */
	memcpy((*cp_data)->x,  x,         dsize);
	memcpy((*cp_data)->xp, sdata->xp, dsize);

	/*store cur_algorithm state*/
	(*cp_data)->algorithm=sdata->cur_algorithm;

	return 0;
}


/*************************************************************************/

int ana_restore_solver_check_point(
		  sca_solv_data* sdata,		            /**< internal solver data */
		  double* x,					        /**< solution vector */
		  sca_solv_checkpoint_data** cp_data    /**< internal check point data */
		  )
{
	unsigned long size;

	if(sdata==NULL) return 1;

	size=sdata->size;

	/*sizes incompatible*/
	if(((*cp_data)==NULL) || (size!=(*cp_data)->size)) return 2;

	/*restore xp and x*/
	memcpy(x,        (*cp_data)->x,  size*sizeof(double));
	memcpy(sdata->xp,(*cp_data)->xp, size*sizeof(double));

	sdata->cur_algorithm=(*cp_data)->algorithm;

	/*printf("Restore solver checkpoint   alg: %i\n",(*cp_data)->cur_algorithm);*/

	return 0;
}

/**************************************************************************/

int ana_copy_solver_check_point(
		  sca_solv_checkpoint_data** cp_data_dest,    /**< internal check point data destination*/
		  sca_solv_checkpoint_data*  cp_data_source   /**< internal check point data source*/
		  )
{
	unsigned long size;

	if(cp_data_source==NULL) return 1;

	size=cp_data_source->size;

	if((*cp_data_dest)==NULL)
	{
		(*cp_data_dest)=(sca_solv_checkpoint_data*)malloc(sizeof(sca_solv_checkpoint_data));
		if((*cp_data_dest)==NULL) return 1;
		(*cp_data_dest)->size=size;
		(*cp_data_dest)->x=(double*)malloc(size*sizeof(double));
		(*cp_data_dest)->xp=(double*)malloc(size*sizeof(double));
		(*cp_data_dest)->algorithm=cp_data_source->algorithm;

		if((*cp_data_dest)->x==NULL) return 1;
		if((*cp_data_dest)->xp==NULL) return 1;
	}

	if(size!=(*cp_data_dest)->size)                 return 2;

	memcpy((*cp_data_dest)->x ,cp_data_source->x, size*sizeof(double));
	memcpy((*cp_data_dest)->xp,cp_data_source->xp,size*sizeof(double));

	(*cp_data_dest)->algorithm=cp_data_source->algorithm;

	return 0;
}

/**************************************************************************/

int ana_allocate_solver_check_point(
		  unsigned long size,                   /**< number of state variables */
		  sca_solv_checkpoint_data** cp_data    /**< internal check point data */
		  )
{
	(*cp_data)=(sca_solv_checkpoint_data*)malloc(sizeof(sca_solv_checkpoint_data));
	if((*cp_data)==NULL) return 1;

	(*cp_data)->x=(double*)calloc(size,sizeof(double));
	if((*cp_data)->x==NULL)
	{
		free(*cp_data);
		(*cp_data)=NULL;
		return 1;
	}

	(*cp_data)->xp=(double*)calloc(size,sizeof(double));
	if((*cp_data)->x==NULL)
	{
		free((*cp_data)->x);
		free(*cp_data);
		(*cp_data)=NULL;
		return 1;
	}

	(*cp_data)->size=size;
	(*cp_data)->algorithm=EULER;

	return 0;
}


/**************************************************************************/

void free_checkpoint_data(sca_solv_checkpoint_data** cp_data)
{
	if((*cp_data)!=NULL)
	{
		if((*cp_data)->x!=NULL)  free((*cp_data)->x);
		(*cp_data)->x=NULL;
		if((*cp_data)->xp!=NULL) free((*cp_data)->xp);
		(*cp_data)->xp=NULL;

		free (*cp_data);
		(*cp_data)=NULL;
	}
}


/****************************************/
