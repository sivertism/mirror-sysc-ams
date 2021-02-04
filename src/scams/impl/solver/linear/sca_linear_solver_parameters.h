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

 sca_linear_solver_parameters.h - define methods set_ignore_woodbury and reset_ignore_woodbury - refer to sca_linear_solver.cpp

 Original Author: Christiane Reuther Fraunhofer IIS/EAS Dresden

 Created on: 27.01.2012

 SVN Version       :  $Revision: 1715 $
 SVN last checkin  :  $Date: 2014-04-27 16:37:25 +0200 (Sun, 27 Apr 2014) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: sca_linear_solver_parameters.h 1715 2014-04-27 14:37:25Z karsten $

 *****************************************************************************/

#ifndef SCA_LINEAR_SOLVER_PARAMETERS_H_
#define SCA_LINEAR_SOLVER_PARAMETERS_H_


namespace sca_core
{
	namespace sca_implementation
	{
		/**
		 * \brief The factorization of coefficient matrices in systems of linear
		 * equations will be done by reinitialization using Gaussian
		 * elimination with privotal search. The cur_algorithm of Woodbury will be
		 * not used.
		 */
		void set_ignore_woodbury();

		/**
		 * \brief The factorization of coefficient matrices in systems of linear
		 * equations will be done by reinitialization using Gaussian
		 * elimination with privotal search and the cur_algorithm of Woodbury.
		 */
		void reset_ignore_woodbury();

		/**
		 * \brief forces the application of the implicit Euler method for
		 * solving the linear DAE (no trapezoidal method)
		 */
		void set_force_implicit_euler_method();

		/**
		 * \brief forces the application of the implicit Euler method and the
		 * trapezoidal method for solving the linear DAE
		 */
		void reset_force_implicit_euler_method();
	}
}

#endif /* SCA_LINEAR_SOLVER_PARAMETERS_H_ */
