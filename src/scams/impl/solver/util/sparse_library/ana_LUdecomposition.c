/*****************************************************************************

    Copyright 2010 - 2011
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

 ana_LUdecomposition.c - description

 Original Author: Christiane Reuther Fraunhofer IIS/EAS Dresden

 Created on: 16.12.2011

 SVN Version       :  $Revision: 1754 $
 SVN last checkin  :  $Date: 2014-05-26 11:18:44 +0200 (Mon, 26 May 2014) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: ana_LUdecomposition.c 1754 2014-05-26 09:18:44Z karsten $

 *****************************************************************************/

/**
 * @file 	ana_LUdecomposition.c
 * @author	Christiane Reuther
 * @date	November 07, 2012
 * @brief	Source-file to define methods <i>ana_LU_decomp_A</i>,
 * <i>ana_get_rank_A</i>, <i>ana_get_size_L</i>, <i>ana_get_size_R</i>
 */

/*****************************************************************************/

#include "ana_solv_data.h"
#include "ma_typedef.h"
#include "ma_util.h"
#include "ma_sparse.h"
#include <string.h>
#include <stdlib.h>

#include "linear_analog_solver.h" /*contains prototype declarations*/

/*****************************************************************************/

/**
 * The method <i>ana_LU_decomp_A</i> generates triangular matrices \f$L\f$ and
 * \f$R\f$ such that \f$A = L\,R\f$ except for permutation. These matrices are
 * stored in <i>(*sdatap)->L</i> and <i>(*sdatap)->R</i>.
 *
 * The factorization of the singular matrix \f$A\f$ is carried out by the
 * function <i>MA_LUdecomposition</i>.
 *
 * In case of <i>first_dec</i> = 1: allocate memory for triangular matrices
 * \f$L\f$ and \f$R\f$.
 *
 * @return
 *  <ul><li>    0 - okay
 *  <li>        1 - reallocation of NULL pointer
 *  <li>        2 - not enough memory
 *  <li>		3 - dimension erroneous
 *  </ul>
 */
int ana_LU_decomp_A_sparse(sca_solv_data **sdatap, int first_dec)
{
	sca_solv_data *sdata;

	struct sparse* sA;
	struct spcode* code;
	/*sparse matrix B equals sparse matrix A, but non-zero entries are ordered*/

	int err = 0;

	sdata = *sdatap;

	if (!first_dec)
	{
		MA_FreeTriangular(sdata->L);
		MA_FreeTriangular(sdata->R);

		free(sdata->L);
		free(sdata->R);
	}

	sA = (struct sparse*) calloc (1, (unsigned) sizeof(struct sparse));
	sdata->L = (struct triangular*) calloc (1,
			(unsigned) sizeof(struct triangular));
	sdata->R = (struct triangular*) calloc (1,
			(unsigned) sizeof(struct triangular));
	code = (struct spcode*) calloc (1, (unsigned) sizeof(struct spcode));

	if (sA == NULL || sdata->L == NULL || sdata->R == NULL || code == NULL)
	{
		if(code!=NULL)      free(code);
		code=NULL;

		if(sdata->R!=NULL)  free(sdata->R);
		sdata->R=NULL;

		if(sdata->L!=NULL)  free(sdata->L);
		sdata->L=NULL;

		if(sA!=NULL) free(sA);
		sA=NULL;

		return 2;
	}

	MA_InitSparse(sA);
	err = MA_GenerateProductValueSparse(sA, sdata->A, 1.0);
	/*sets sparse matrix entries in correct order for factorization */
	if (err)	goto end;

	MA_InitCode(code);
	MA_InitTriangular(sdata->L);
	MA_InitTriangular(sdata->R);

	sdata->L->a = (value *) calloc((unsigned) sA->nel,
			(unsigned) sizeof(value));
	sdata->L->column = (count_near *) calloc((unsigned) sA->nel,
			(unsigned) sizeof(count_near));
	sdata->L->line = (count_near *) calloc((unsigned) sA->nel,
			(unsigned) sizeof(count_near));

	sdata->R->a = (value *) calloc((unsigned) sA->nmax,
			(unsigned) sizeof(value));
	sdata->R->column = (count_near *) calloc((unsigned) sA->nmax,
			(unsigned) sizeof(count_near));
	sdata->R->line = (count_near *) calloc((unsigned) sA->nmax,
			(unsigned) sizeof(count_near));

	if (sdata->L->a == NULL || sdata->L->column == NULL
			|| sdata->L->line == NULL || sdata->R->a == NULL
			|| sdata->R->column == NULL || sdata->R->line == NULL)
	{
		err = 2;
		goto end;
	}

	err = MA_LUdecomposition(sA, code, sdata->L, sdata->R);
	if (err == 0)
		sdata->rank_A = code->rank;

	end:

	MA_FreeSparse(sA);
	MA_FreeCode(code);

	free(sA);
	free(code);

	return err;
}

/*****************************************************************************/

/**
 * @return
 *  <ul><li>    0 - okay
 *  <li>        1 - reallocation of NULL pointer
 *  <li>        2 - not enough memory
 *  <li>		3 - dimension erroneous
 *  </ul>
 */
int ana_LU_decomp_A(sca_solv_data **sdatap, int first_dec)
{
	int err = 0;

	err = ana_LU_decomp_A_sparse(sdatap, first_dec);

	return err;
}

/*****************************************************************************/

/**
 * The method <i>ana_get_rank_A</i> outputs the rank of the matrix \f$A\f$
 * <i>data->rank_A</i>.
 */
long ana_get_rank_A(sca_solv_data* data)
{
	return data->rank_A;
}

/*****************************************************************************/

/**
 * The method <i>ana_get_size_L</i> outputs the number of non-zero elements of
 * the triangular matrix \f$L\f$.
 */
long ana_get_size_L(sca_solv_data* data)
{
	return data->L->nel;
}

/*****************************************************************************/

/**
 * The method <i>ana_get_size_R</i> outputs the number of non-zero elements of
 * the triangular matrix \f$R\f$.
 */
long ana_get_size_R(sca_solv_data* data)
{
	return data->R->nel;
}

/*****************************************************************************/
