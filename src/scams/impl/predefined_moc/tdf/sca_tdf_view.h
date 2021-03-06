/*****************************************************************************

    Copyright 2010
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

 sca_tdf_view.h - description

 Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

 Created on: 05.08.2009

 SVN Version       :  $Revision: 1265 $
 SVN last checkin  :  $Date: 2011-11-19 20:43:31 +0000 (Sat, 19 Nov 2011) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: sca_tdf_view.h 1265 2011-11-19 20:43:31Z karsten $

 *****************************************************************************/

/*****************************************************************************/

#ifndef SCA_TDF_VIEW_H_
#define SCA_TDF_VIEW_H_

#include "scams/impl/core/sca_view_base.h"

namespace sca_tdf
{

namespace sca_implementation
{

/**
 timed static dataflow view implementation
 */

class sca_tdf_view: public sca_core::sca_implementation::sca_view_base
{

public:

	/** Overwritten (virtual) method for equation setup */
	void setup_equations();

	sca_tdf_view();
	~sca_tdf_view();

};

} //namespace sca_implementation
}//namespace sca_tdf


#endif /* SCA_TDF_VIEW_H_ */
