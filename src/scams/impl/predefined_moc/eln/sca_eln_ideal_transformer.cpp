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

  sca_eln_ideal_transformer.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

  Created on: 10.11.2009

   SVN Version       :  $Revision: 1910 $
   SVN last checkin  :  $Date: 2016-02-16 13:44:37 +0100 (Tue, 16 Feb 2016) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_ideal_transformer.cpp 1910 2016-02-16 12:44:37Z karsten $

 *****************************************************************************/

/*****************************************************************************/


#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_ideal_transformer.h"

namespace sca_eln
{

sca_ideal_transformer::sca_ideal_transformer(sc_core::sc_module_name, double ratio_) :
	p1("p1"), n1("n1"), p2("p2"), n2("n2"), ratio("ratio", ratio_)
{
	through_value_available = false;
	through_value_type = "I";
	through_value_unit = "A";

	nadd=-1;
}


const char* sca_ideal_transformer::kind() const
{
	return "sca_eln::sca_ideal_transformer";
}

void sca_ideal_transformer::matrix_stamps()
{
    nadd = add_equation();

    B_wr(nadd, p1) =  1.0;
    B_wr(nadd, n1) = -1.0;
    B_wr(nadd, p2) = -ratio;
    B_wr(nadd, n2) =  ratio;
    B_wr(p1,   nadd) =  1.0;
    B_wr(n1,   nadd) = -1.0;
    B_wr(p2,   nadd) = -ratio;
    B_wr(n2,   nadd) =  ratio;
}

} //namespace sca_eln

