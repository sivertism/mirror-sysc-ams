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

  sca_eln_tdf_l.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden / COSEDA Technologies

  Created on: 10.11.2009

   SVN Version       :  $Revision: 1892 $
   SVN last checkin  :  $Date: 2016-01-10 12:59:12 +0100 (Sun, 10 Jan 2016) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_tdf_l.cpp 1892 2016-01-10 11:59:12Z karsten $

 *****************************************************************************/

/*****************************************************************************/


#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_tdf_l.h"
#include "scams/impl/util/data_types/sca_method_object.h"
#include "scams/impl/util/data_types/sca_method_list.h"
#include "scams/impl/core/sca_solver_base.h"
#include <limits>
#include "scams/impl/util/data_types/almost_equal.h"
#include "scams/impl/predefined_moc/conservative/sca_con_interactive_trace_data.h"


namespace sca_eln
{

namespace sca_tdf
{

sca_l::sca_l(sc_core::sc_module_name, double scale_, double psi0_) :
p("p"), n("n"), inp("inp"), scale("scale", scale_), psi0("psi0",psi0_), phi0(psi0)
{
    through_value_available = true;
    through_value_type      = "I";
    through_value_unit      = "A";

    curr_value=1.0;

    nadd1=-1;
    nadd2= -1;

    dc_init=true;
}

const char* sca_l::kind() const
{
	return "sca_eln::sca_tdf::sca_l";
}

void sca_l::trace( sc_core::sc_trace_file* tf ) const
{
	std::ostringstream str;
	str << "sc_trace of sca_ln::sca_tdf::sca_l module not supported for module: " << this->name();
	SC_REPORT_WARNING("SystemC-AMS",str.str().c_str());
}

void sca_l::read_lval()
{

    double inp_val = scale * inp.read();

	if(!sca_util::sca_implementation::almost_equal(inp_val,curr_value))
    {
    	enable_b_change=true;
    	continous=true;

		if (!sca_ac_analysis::sca_ac_is_running() && dc_init &&
	       			(psi0.get()== sca_util::SCA_UNDEFINED))
		{
			B_wr(nadd2, nadd1) = inp_val;
		}
		else
		{
			B_wr(nadd1, nadd1) = inp_val;
		}

		curr_value=inp_val;

    }
}



void sca_l::matrix_stamps()
{
	// nadd1 - i
	// nadd2 - psi
    nadd1 = add_equation();
    nadd2 = add_equation();

    add_method(PRE_SOLVE, SCA_VMPTR(sca_l::read_lval));

	if(!sca_ac_analysis::sca_ac_is_running() && dc_init)
	{
		if (psi0.get() == sca_util::SCA_UNDEFINED)
		{
			//short cut
			B_wr(nadd1, p) = 1.0;
			B_wr(nadd1, n) = -1.0;

			B_wr(p, nadd1) = 1.0;
			B_wr(n, nadd1) = -1.0;

			B_wr(nadd2, nadd2) = -1.0; //psi=L*i
			B_wr(nadd2, nadd1) = curr_value;
		}
		else
		{
			B_wr(p, nadd1) = 1.0;
			B_wr(n, nadd1) = -1.0;

			// L*i = psi
			B_wr(nadd1, nadd1) = curr_value;
			B_wr(nadd1, nadd2) = -1.0;

			// psi = psi0
			B_wr(nadd2, nadd2) = -1.0;
			q(nadd2).set_value(psi0.get());
		}
		add_method(POST_SOLVE, SCA_VMPTR(sca_l::post_solve));
	}
	else
	{
		B_wr(p, nadd1) =  1.0;
		B_wr(n, nadd1) = -1.0;

		// L*i = psi
		B_wr(nadd1, nadd1) =  curr_value;
		B_wr(nadd1, nadd2) = -1.0;

	    //dpsi = v
	    B_wr(nadd2, p) =  1.0;
	    B_wr(nadd2, n) = -1.0;
	    A(nadd2, nadd2) = -1.0;
	}
}


void sca_l::post_solve() //set capacitor stamps after first step
{
	//reset dc stamps
	B_wr(nadd2, nadd2) =  0.0;
	q(nadd2).set_value(0.0);

	B_wr(nadd2, nadd1) = 0.0;

	B_wr(nadd1, p) = 0.0;
	B_wr(nadd1, n) = 0.0;

	B_wr(p, nadd1) = 0.0;
	B_wr(n, nadd1) = 0.0;

	remove_method(POST_SOLVE, SCA_VMPTR(sca_l::post_solve));

	//set transient stamps
	B_wr(p, nadd1) =  1.0;
	B_wr(n, nadd1) = -1.0;

	// L*i = psi
	B_wr(nadd1, nadd1) =  curr_value;
	B_wr(nadd1, nadd2) = -1.0;

    //dpsi = v
    B_wr(nadd2, p) =  1.0;
    B_wr(nadd2, n) = -1.0;
    A(nadd2, nadd2) = -1.0;

    request_reinit(1);

	dc_init=false;

}



bool sca_l::trace_init(sca_util::sca_implementation::sca_trace_object_data& data)
{
    data.type=through_value_type;
    data.unit=through_value_unit;

    //trace will be activated after every complete cluster calculation
    //by teh synchronization layer
    get_sync_domain()->add_solver_trace(data);
    return true;
}

void sca_l::trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)
{
    sca_core::sca_time ctime = sca_eln::sca_module::get_time();

    double through_value = x(nadd1);

	if((id<0)&&(this->trd!=NULL)) //interactive trace
	{
	    this->trd->store_value(through_value);
		return;
	}

    buffer.store_time_stamp(id,ctime,through_value);
}

sca_util::sca_complex sca_l::calculate_ac_result(sca_util::sca_complex* res_vec)
{
	return res_vec[nadd1];
}



const double& sca_l::get_typed_trace_value() const
{
	if(this->trd==NULL)
	{
		this->trd=new sca_core::sca_implementation::sca_con_interactive_trace_data(this);

	}

	return this->trd->get_value();
}


const std::string& sca_l::get_trace_value() const
{
	std::ostringstream str;
	str << this->get_typed_trace_value();

	static std::string trace_value_string;
	trace_value_string=str.str().c_str();
	return trace_value_string;
}





bool sca_l::register_trace_callback(sca_trace_callback cb,void* cb_arg)
{
	if(this->trd==NULL)
	{
		this->trd=new sca_core::sca_implementation::sca_con_interactive_trace_data(this);
	}


	this->trd->register_trace_callback(cb,cb_arg);

	return true;

}


} //namespace sca_tdf
} //namespace sca_eln

