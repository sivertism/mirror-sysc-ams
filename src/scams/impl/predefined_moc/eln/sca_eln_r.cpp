/*****************************************************************************

    Copyright 2010-2013
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

  sca_eln_r.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden / COSEDA Technologies

  Created on: 10.11.2009

   SVN Version       :  $Revision: 1892 $
   SVN last checkin  :  $Date: 2016-01-10 12:59:12 +0100 (Sun, 10 Jan 2016) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_r.cpp 1892 2016-01-10 11:59:12Z karsten $

 *****************************************************************************/

/*****************************************************************************/

#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_r.h"
#include "scams/impl/core/sca_solver_base.h"
#include "scams/impl/predefined_moc/conservative/sca_con_interactive_trace_data.h"

namespace sca_eln
{
sca_r::sca_r(sc_core::sc_module_name, double value_) :
		p("p"), n("n"), value("value", value_)
{
    through_value_available = true;
    through_value_type      = "Current";
    through_value_unit      = "A";
    nadd=-1;
}

const char* sca_r::kind() const
{
	return "sca_eln::sca_r";
}


void sca_r::trace( sc_core::sc_trace_file* tf ) const
{
	std::ostringstream str;
	str << "sc_trace of sca_r module not supported for module: " << this->name();
	SC_REPORT_WARNING("SystemC-AMS",str.str().c_str());
}


void sca_r::matrix_stamps()
{
	if(fabs(value.get())>1e-4)
	{
		B_wr(p,p) +=  1.0/value.get();
		B_wr(p,n) += -1.0/value.get();
		B_wr(n,p) += -1.0/value.get();
		B_wr(n,n) +=  1.0/value.get();
	}
	else
	{
	    nadd = add_equation();

	    B_wr(nadd, p) =  1.0;
	    B_wr(nadd, n) = -1.0;
	    B_wr(p, nadd) =  1.0;
	    B_wr(n, nadd) = -1.0;

	    B_wr(nadd, nadd) = -value.get();
	}
}

bool sca_r::trace_init(sca_util::sca_implementation::sca_trace_object_data& data)
{
    data.type=through_value_type;
    data.unit=through_value_unit;

    //trace will be activated after every complete cluster calculation
    //by teh synchronization layer
    get_sync_domain()->add_solver_trace(data);
    return true;
}

void sca_r::trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)
{
    sca_core::sca_time ctime = sca_eln::sca_module::get_time();

    double through_value;

    if(nadd<0)
    {
    	through_value = (x(p)-x(n))/value.get();
    }
    else
    {
    	through_value=x(nadd);
    }


	if((id<0)&&(this->trd!=NULL)) //interactive trace
	{
	    this->trd->store_value(through_value);
		return;
	}

    buffer.store_time_stamp(id,ctime,through_value);
}

sca_util::sca_complex sca_r::calculate_ac_result(sca_util::sca_complex* res_vec)
{
    //if reference node return 0.0
    sca_util::sca_complex rp = long(p) < 0 ? 0.0 : res_vec[p];
    sca_util::sca_complex rn = long(n) < 0 ? 0.0 : res_vec[n];

    if(nadd>=0) return res_vec[nadd];

    return (rp-rn)/value.get();
}



const double& sca_r::get_typed_trace_value() const
{
	if(this->trd==NULL)
	{
		this->trd=new sca_core::sca_implementation::sca_con_interactive_trace_data(this);

	}

	return this->trd->get_value();
}


const std::string& sca_r::get_trace_value() const
{
	std::ostringstream str;
	str << this->get_typed_trace_value();

	static std::string trace_value_string;
	trace_value_string=str.str().c_str();
	return trace_value_string;
}





bool sca_r::register_trace_callback(sca_trace_callback cb,void* cb_arg)
{
	if(this->trd==NULL)
	{
		this->trd=new sca_core::sca_implementation::sca_con_interactive_trace_data(this);
	}


	this->trd->register_trace_callback(cb,cb_arg);

	return true;

}



} //namespace sca_eln
