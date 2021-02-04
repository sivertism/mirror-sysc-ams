/*****************************************************************************

    Copyright 2010-2013
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

  sca_eln_tdf_isink.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

  Created on: 10.11.2009

   SVN Version       :  $Revision: 1910 $
   SVN last checkin  :  $Date: 2016-02-16 13:44:37 +0100 (Tue, 16 Feb 2016) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_tdf_isink.cpp 1910 2016-02-16 12:44:37Z karsten $

 *****************************************************************************/

/*****************************************************************************/


#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_tdf_isink.h"
#include "scams/impl/util/data_types/sca_method_object.h"
#include "scams/impl/core/sca_solver_base.h"

namespace sca_eln
{

namespace sca_tdf
{


const char* sca_isink::kind() const
{
	return "sca_eln::sca_tdf::sca_isink";
}


void sca_isink::trace( sc_core::sc_trace_file* tf ) const
{
	std::ostringstream str;
	str << "sc_trace of sca_ln::sca_tdf::sca_isink module not supported for module: " << this->name();
	SC_REPORT_WARNING("SystemC-AMS",str.str().c_str());
}


sca_isink::	sca_isink(sc_core::sc_module_name, double scale_) :
	p("p"), n("n"), outp("outp"), scale("scale", scale_)
{
    through_value_available = true;
    through_value_type      = "I";
    through_value_unit      = "A";

    nadd=-1;
}


void sca_isink::assign_result()
{

    if(!sca_ac_analysis::sca_ac_is_running())
    {
    	outp.write(x(nadd) * scale );
    }
    else
    {
        sca_ac_analysis::sca_ac(outp) = scale * x(nadd);
    }
}

void sca_isink::matrix_stamps()
{
    nadd = add_equation();

    B_wr(nadd,p) +=  1.0;
    B_wr(nadd,n) += -1.0;
    B_wr(p,nadd) +=  1.0;
    B_wr(n,nadd) += -1.0;

    add_method(POST_SOLVE, SCA_VMPTR(sca_isink::assign_result));
}

bool sca_isink::trace_init(sca_util::sca_implementation::sca_trace_object_data& data)
{
    data.type=through_value_type;
    data.unit=through_value_unit;

    //trace will be activated after every complete cluster calculation
    //by teh synchronization layer
    get_sync_domain()->add_solver_trace(data);
    return true;
}

void sca_isink::trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)
{
    sca_core::sca_time ctime = sca_eln::sca_module::get_time();

    double through_value = x(nadd);
    buffer.store_time_stamp(id,ctime,through_value);
}

sca_util::sca_complex sca_isink::calculate_ac_result(sca_util::sca_complex* res_vec)
{
	return res_vec[nadd];
}



} //namespace sca_tdf
} //namespace sca_eln

