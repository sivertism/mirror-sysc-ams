/*****************************************************************************

    Copyright 2012
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

 ma_woodbury.h - description

 Original Author: Christiane Reuther, Fraunhofer IIS/EAS Dresden

 Created on: September 24, 2012

 SVN Version       :  $Revision: 1401 $
 SVN last checkin  :  $Date: 2012-11-20 14:38:10 +0100 (Tue, 20 Nov 2012) $
 SVN checkin by    :  $Author: reuther $
 SVN Id            :  $Id: ma_woodbury.h 1401 2012-11-20 13:38:10Z reuther $

 *****************************************************************************/

/**
 * @file 	ma_woodbury.h
 * @author	Christiane Reuther
 * @date	September 24, 2012
 * @brief	Header-file to define method for the formula of Woodbury
 */

/*****************************************************************************/

#ifndef _ma_woodbury_h_
#define _ma_woodbury_h_

/* /// Headers ////////////////////////////////////////////////////////////// */

#include "ma_typedef.h"
#include "ana_solv_data.h"

/* /// Defines ////////////////////////////////////////////////////////////// */

#define exportMA_Sparse

/*MA_woodbury.c*/

/**
 * \brief computes solution of linear system of equations by using the
 * formula of Woodbury
 */
exportMA_Sparse err_code MA_LowRankModifications(
		sca_solv_data* sdata,			/**< internal solver data */
		value* x,						/**< solution vector */
		count_far k,					/**< number of modifications */
		count_near* mode,				/**< kinds of modification */
		count_far* n,					/**< positions of modification */
		value* s						/**< values of modification */
		);

#endif
