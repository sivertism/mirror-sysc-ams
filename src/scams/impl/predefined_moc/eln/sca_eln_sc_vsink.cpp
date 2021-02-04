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

  sca_eln_sc_vsink.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

  Created on: 10.11.2009

   SVN Version       :  $Revision: 1665 $
   SVN last checkin  :  $Date: 2013-12-23 23:12:50 +0100 (Mon, 23 Dec 2013) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_sc_vsink.cpp 1665 2013-12-23 22:12:50Z karsten $

 *****************************************************************************/

/*****************************************************************************/


#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_sc_vsink.h"
#include "scams/impl/util/data_types/sca_method_object.h"
#include "scams/impl/core/sca_solver_base.h"

namespace sca_eln
{

namespace sca_de
{


const char* sca_vsink::kind() const
{
	return "sca_eln::sca_de::sca_vsink";
}


void sca_vsink::trace( sc_core::sc_trace_file* tf ) const
{
	std::ostringstream str;
	str << "sc_trace of sca_ln::sca_de::sca_vsink module not supported for module: " << this->name();
	SC_REPORT_WARNING("SystemC-AMS",str.str().c_str());
}




sca_vsink::sca_vsink(sc_core::sc_module_name, double scale_) :
	p("p"), n("n"), outp("outp"), scale("scale", scale_)
{

    //IMPROVE: find better solution -> inserts additional port in database
    conv_port=new ::sca_tdf::sca_de::sca_out<double>("converter_port");
    conv_port->bind(outp);
}


void sca_vsink::assign_result()
{

    if(!sca_ac_analysis::sca_ac_is_running())
    {
    	conv_port->write((x(p)-x(n))*scale);
    }
    else
    {
        sca_ac_analysis::sca_implementation::sca_ac(outp) = scale * x(p);
    }

}

void sca_vsink::matrix_stamps()
{
	add_method(POST_SOLVE, SCA_VMPTR(sca_vsink::assign_result));
}

bool sca_vsink::trace_init(sca_util::sca_implementation::sca_trace_object_data& data)
{
    data.type=through_value_type;
    data.unit=through_value_unit;

    //trace will be activated after every complete cluster calculation
    //by teh synchronization layer
    get_sync_domain()->add_solver_trace(data);
    return true;
}
void sca_vsink::trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)
{
    sca_core::sca_time ctime = sca_eln::sca_module::get_time();

    double through_value = 0.0;
    buffer.store_time_stamp(id,ctime,through_value);
}

sca_util::sca_complex sca_vsink::calculate_ac_result(sca_util::sca_complex* res_vec)
{
	return sca_util::sca_complex(0.0,0.0);
}



} //namespace sca_de
} //namespace sca_eln
