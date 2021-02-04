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

 MA_woodbury.c - description

 Original Author: Christiane Reuther Fraunhofer IIS/EAS Dresden

 Created on:

 SVN Version       :  $Revision: 1754 $
 SVN last checkin  :  $Date: 2014-05-26 11:18:44 +0200 (Mon, 26 May 2014) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: MA_woodbury.c 1754 2014-05-26 09:18:44Z karsten $

 *****************************************************************************/

/**
 * @file 	MA_woodbury.c
 * @author	Christiane Reuther
 * @date	September 24, 2012
 * @brief	Source-file to define method for the formula of Woodbury
 */

/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ma_woodbury.h"
#include "ma_sparse.h"

/*////////////////////////////////////////////////////////////////////////*/

/**
 *  Detailed information about the formula of Woodbury within the linear solver
 *  of SystemC AMS is given in the paper [Reuther, Einwich: A SystemC AMS
 *  extension for controlled modules and dynamic step sizes, Forum on
 *  Specification and Design Languages (FDL) 2012 , 2012, 90 - 97].
 *
 *	The formula of Woodbury is as follows: <br>
 *	Let the matrices \f$ A \f$ and \f$ I + W A^{-1} V \f$ be square nonsingular
 *	matrices and the matrices \f$ V \f$ and \f$ W \f$ rectangular and \f$ I \f$
 *	the unit matrix such that all appearing matrix multiplications are valid.
 *	Then the matrix \f$ A + V W \f$ is nonsingular and \f$(A + V W)^{-1} =
 *	A^{-1} - A^{-1} V \left(I + W A^{-1} V\right)^{-1} W A^{-1} \f$.
 *
 *	Implementation details of the function MA_LowRankModifications: <br>
 *	<ul>
 *	<li> The vector \f$ x \f$ deals as input by \f$ x = A^{-1} r \f$, where
 *	\f$ r \f$ is the right-hand-side of the equation system to be solved.
 *	<li> The vector \f$ x \f$ deals as output by \f$ x = (A + V W)^{-1} r \f$
 *	as the solution of the equation system.
 *	<li> The matrix \f$ A \f$ is given by <i>sdata->sZ_euler</i> or
 *	<i>sdata->sZ_trapez</i> in dependence on the choice of the solution
 *	cur_algorithm <i>sdata->cur_algorithm</i> in the linear solver.
 *	<li> In case of <i>k</i> > 0 the matrices \f$ V \f$ and \f$ W \f$ are given
 *	by the parameters <i>mode, n, s</i>, where
 *	<ul>
 *	<li> <i>mode</i> = 1: 1 perturbation on diagonal in \f$ A \f$ - unit vectors
 *	 in	\f$ V \f$ and \f$ W \f$
 *	<li> <i>mode</i> = 2: 1 perturbation not on diagonal in \f$ A \f$ - unit
 *	vectors in \f$ V \f$ and \f$ W \f$
 *	<li> <i>mode</i> = 3: 2 perturbations (same absolute value) in \f$ A \f$ in
 *	one line - sum of two unit vectors in \f$ V \f$, unit vector in \f$ W \f$
 *	<li> <i>mode</i> = 4: 4 perturbations (same absolute value) in \f$ A \f$ in
 *	two lines and two colomns - sum of two unit vectors in \f$ V \f$ and
 *	\f$ W \f$
 *	<li> <i>mode</i> = 5: 2 perturbations (same absolute value) in \f$ A \f$ in
 *	one column - sum of two unit vectors in \f$ W \f$, unit vector in \f$ V \f$
 *	</ul>
 *	and <i>s</i> is the vector of values of the modifications and <i>n</i>
 *	contains the positions of the modifications.
 *	<li> In case of <i>sdata->variable_step_size</i> = 1 and the rank of the
 *	matrix <i>sdata->sW</i> > 0 the matrix \f$ V \f$ contains the triangular
 *	matrix <i>sdata->L</i> and the matrix \f$ W \f$ contains the triangular
 *	matrix <i>sdata->R</i>, where \f$ W = L R \f$ without mentioning necessary
 *	permutation matrices.
 *	<li> The matrix \f$V\f$ is of size \f$dim_A \times rank_{wb}\f$ and the
 *	matrix \f$W\f$ is of size \f$rank_{wb} \times dim_A\f$, where \f$dim_A\f$ is
 *	the dimension of \f$A\f$ and \f$rank_{wb}\f$ is the sum of \f$k\f$ and the
 *	rank of	the matrix \f$W\f$.
 *	</ul>
 *
 *  @return
 *  <ul><li>    0 - okay
 *  <li>        1 - reallocation of NULL pointer
 *  <li>        2 - not enough memory
 *  <li>		3 - dimension erroneous
 *  <li>		4 - matrix singular
 *  <li>		6 - no sparse matrix
 *  </ul>
 */

exportMA_Sparse err_code MA_LowRankModifications(
		sca_solv_data* sdata,
		value* x,
		count_far k,
		count_near* mode,
		count_far* n,
		value* s
		)
{
	/* mode == 1: 1 perturbation on diagonal - 1 entry in n = position on
	 * diagonal - e.g. sca_eln_tdf_r.cpp */
	/* mode == 2: 1 perturbation not on diagonal - 2 entries in n = position
	 * in matrix - e.g. sca_eln_tdf_l.cpp */
	/*  mode == 3: 2 perturbations (same absolute value_t) in matrix in one line
	 * - 3 entries in n = position of line and positions of columns
	 * - e.g. sca_eln_tdf_c.cpp */
	/* mode == 4: 4 perturbations (same absolute value_t) in matrix in two lines
	 * and two colomns - 2 entries in n = positions of lines and positions of
	 * columns, resp. - e.g. sca_eln_tdf_rswitch.cpp */
	/* mode == 5: 2 perturbations (same absolute value_t) in matrix in one
	 * column - 3 entries in n = positions of lines and position of column
	 * - e.g. sca_lsf_tdf_demux.cpp */

	value* v; /* vector with dim entries - column vectors of matrix V*/
	value* xv; /* matrix with dim*rank_wb entries - xv = A^{-1}*v */
	value* v1; /* matrix with rank_wb*rank_wb entries - v1 = I + W*xv */
	value* v3; /* vector with rank_wb entries - v3 = W*x */
	value* v5; /* vector with rank_wb entries - v5 = v1^{-1}*v3 */
	value* v2; /* vector with dim entries - v2 = xv*v5 */

	struct sparse* ksparse=NULL;
	struct spcode* kcode=NULL;

	struct sparse* sZ = NULL;
	struct spcode* code = NULL;

	count_far dim, r;

	int i0, j0, i, j, m;
	int err = 0;
	long rank_wb;

    if (sdata->cur_algorithm == EULER)
    {
    	sZ = sdata->sZ_euler;
    	code = sdata->code_euler;
    }
    else /*if (sdata->cur_algorithm == TRAPEZ)*/
    {
    	sZ = sdata->sZ_trapez;
    	code = sdata->code_trapez;
    }

    if ((sdata->variable_step_size == 0) || (sdata->rank_A == 0))
    	r = 0;
	else
    	r = sdata->rank_A;

	dim = sZ->nd;

	if (fabs(sdata->h_diff) < 1.0e-18)
		r = 0;

	rank_wb = r + k;

	if (rank_wb == 0)
		return 0;

	memset(sdata->aux_woodbury, 0, 2*sZ->nd*sZ->nd * sizeof(double));

	v = sdata->aux_woodbury;
	xv = v + dim;
	v1 = xv + dim*rank_wb;
	v2 = v1 + rank_wb*rank_wb;
	v3 = v2 + dim;
	v5 = v3 + rank_wb;

/* Woodbury cur_algorithm *******************  */

	/* modification in matrix B */

	for (i = 0; i < k; i++)
	{
		if (mode[i] == 1 || mode[i] == 2 || mode[i] == 3)
		{
			*(v+n[3*i]) = *(s+i); /*n[i]th unit vector*/
		}
		else if (mode[i] == 4 || mode[i] == 5) /*v = ...*/
		{
			*(v + n[3*i]) = *(s+i);
			*(v + n[3*i+1]) = -*(s+i);
		}

		MA_LequSparseSolut(sZ, code, v, v2); /*A*xv = v */

		for (j = 0; j< dim; j++)
		{
			*(xv + i*dim + j) = *(v2 + j);

			*(v+j) = 0;
			*(v2+j) = 0;
		}
	}

	/* variable time step */

	m = 0;

	for (i = 0; i < r; i++)
	{
		while(i == sdata->L->column[m] && m < sdata->L->nel)
		{
			*(v+sdata->L->line[m]) = sdata->h_diff * sdata->L->a[m];
			m++;
		}

		MA_LequSparseSolut(sZ, code, v, v2); /*A*xv = v */

		for (j = 0; j< dim; j++)
		{
			*(xv + (i+k)*dim + j) = *(v2 + j);

			*(v+j) = 0;
			*(v2+j) = 0;
		}
	}

	/* modification in matrix B */

	for (i = 0; i< k; i++)
	{
		if (mode[i] == 1)
		{
			for (j = 0; j < rank_wb; j++) /*v1 = w^T * xv */
			{
				*(v1 +j*rank_wb + i) = *(xv + n[3*i] + j* dim);
			}

			*(v3+i) = *(x+n[3*i]); /*v3 = w^T * x */
		}

		if (mode[i] == 2)
		{
			for (j = 0; j < rank_wb; j++) /*v1 = w^T * xv */
			{
				*(v1 +j*rank_wb + i) = *(xv + n[3*i+1] + j* dim);
			}

			*(v3+i) = *(x+n[3*i+1]); /*v3 = w^T * x */
		}

		if (mode[i] == 3)
		{
			for (j = 0; j < rank_wb; j++) /*v1 = w^T * xv */
			{
				*(v1 +j*rank_wb + i) =
						*(xv + n[3*i+1] + j* dim) - *(xv + n[3*i+2] + j* dim);
			}

			*(v3+i) = *(x+n[3*i+1]) - *(x+n[3*i+2]); /*v3 = w^T * x */
		}

		if (mode[i]== 4)
		{
			for (j = 0; j < rank_wb; j++) /*v1 = w^T * xv */
			{
				*(v1 +j*rank_wb + i) =
						*(xv + n[3*i] + j* dim) - *(xv + n[3*i+1] + j* dim);
			}

			*(v3+i) = *(x+n[3*i]) - *(x+n[3*i+1]); /*v3 = w^T * x */
		}

		if (mode[i] == 5)
		{
			for (j = 0; j < rank_wb; j++) /*v1 = w^T * xv */
			{
				*(v1 +j*rank_wb + i) = *(xv + n[3*i+2] + j* dim);
			}

			*(v3+i) = *(x+n[3*i+2]); /*v3 = w^T * x */
		}

		*(v1 + i*rank_wb + i) += 1; /*v1 = v1 + I */
	}

	/* variable time step */

	m = 0;

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < rank_wb; j++)
		{
			*(v1 + j*rank_wb + k + i) = *(xv + sdata->R->column[m] + j*dim);
		}

		*(v3 + k + i) = *(x+sdata->R->column[m]);
		m++;

		while(i == sdata->R->line[m] && m < sdata->R->nel)
		{
			for (j = 0; j < rank_wb; j++)
			{
				*(v1 + j*rank_wb + k + i) +=
						*(xv + sdata->R->column[m] + j*dim) * sdata->R->a[m];
			}

			*(v3 + k + i) += *(x+sdata->R->column[m]) * sdata->R->a[m];
			m++;
		}

		*(v1 + (k+i)*rank_wb + k + i) += 1; /*v1 = v1 + I */
	}

	/* linear equation system with rank_wb equations */

	if (rank_wb == 1) /* v1*v5 = v3 */
	{
		if (*v1 > code->piv_abs_tol || *v1 < - code->piv_abs_tol)
		{	/*modified matrix is regular */
			*v5 = *v3 / *v1 ;
		}
		else /*modified matrix is singular */
		{
			if (k == 1)
			{
				code->critical_line = *n;

				if (*mode == 1)
					code->critical_column = *n;
				else if (*mode == 5)
					code->critical_column = *(n+2);
				else /*if (*mode >= 1 && *mode <= 5)*/
					code->critical_column = *(n+1);
			}
			else
			{
				code->critical_column = sdata->L->line[0];
				code->critical_line = sdata->R->column[0];
			}

			err = 4;

			goto retour;
		}
	}
	else	/* rank_wb > 1 */
	{
		ksparse = (struct sparse*) calloc(1,(unsigned) sizeof(struct sparse));
		kcode = (struct spcode*) calloc(1,(unsigned) sizeof(struct spcode));

		if (ksparse == NULL || kcode == NULL)
		{
			err = 2;
			goto retour;
		}

		MA_InitSparse(ksparse);
		MA_InitCode(kcode);
		MA_ConvertFullToSparse(v1,rank_wb,ksparse,0);

		kcode->piv_abs_tol = 1.0e-16;

		err = MA_LequSparseCodegen(ksparse, kcode);
		if (err) /*modified matrix is singular*/
		{
			if (!(kcode->critical_column > k))
			{
				code->critical_column = n[3*kcode->critical_column];
			}
			else /*(kcode->critical_column > k)*/
			{
				i0 = 0;
				j0 = sdata->L->column[i0];

				while (j0 < kcode->critical_column - k)
				{
					i0++;
					j0 = sdata->L->column[i0];
				}

				code->critical_column = sdata->L->line[i0];
			}

			if (!(kcode->critical_line > k))
			{
				if (mode[kcode->critical_line] == 1)
					code->critical_line = n[3*kcode->critical_line];
				else if (mode[kcode->critical_line] == 5)
					code->critical_line = n[3*kcode->critical_line+2];
				else /*if (mode[kcode->critical_line] >= 1 && <= 5)*/
					code->critical_line = n[3*kcode->critical_line+1];
			}
			else /*(kcode->critical_line > k)*/
			{
				i0 = 0;
				j0 = sdata->R->line[i0];

				while (j0 < kcode->critical_line - k)
				{
					i0++;
					j0 = sdata->R->line[i0];
				}

				code->critical_line = sdata->R->column[i0];
			}

			goto retour;
		}

		MA_LequSparseSolut(ksparse, kcode, v3, v5);
	}

	for (j = 0; j < dim; j++) /*v2 = xv * v5, x = x - v2 */
	{
		for (i = 0; i < rank_wb; i++)
		{
			*(v2 + j) +=  *(xv + j + i * dim) * *(v5 + i);
		}

		*(x+j) -= *(v2 + j);
	}

/*  Retour *******************  */

	retour:

	if (rank_wb != 1)
	{
		MA_FreeSparse(ksparse);
		free(ksparse);
		MA_FreeCode(kcode);
		free(kcode);
	}

	return err;
}

/* /// end of file ////////////////////////////////////////////////////////// */
