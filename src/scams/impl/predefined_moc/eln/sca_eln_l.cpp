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

  sca_eln_l.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden / COSEDA Technologies

  Created on: 10.11.2009

   SVN Version       :  $Revision: 1863 $
   SVN last checkin  :  $Date: 2016-01-05 10:11:33 +0100 (Tue, 05 Jan 2016) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_l.cpp 1863 2016-01-05 09:11:33Z karsten $

 *****************************************************************************/

/*****************************************************************************/

#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_l.h"
#include "scams/impl/util/data_types/sca_method_object.h"
#include "scams/impl/util/data_types/sca_method_list.h"
#include "scams/impl/core/sca_solver_base.h"
#include "scams/impl/predefined_moc/conservative/sca_con_interactive_trace_data.h"


namespace sca_eln
{

sca_l::sca_l(sc_core::sc_module_name, double value_, double psi0_) :
p("p"), n("n"), value("value", value_), psi0("psi0", psi0_), phi0(psi0)
{
    through_value_available = true;
    through_value_type      = "I";
    through_value_unit      = "A";

    dc_init=true;

    nadd=-1;
}


const char* sca_l::kind() const
{
	return "sca_eln::sca_l";
}

void sca_l::trace( sc_core::sc_trace_file* tf ) const
{
	std::ostringstream str;
	str << "sc_trace of sca_l module not supported for module: " << this->name();
	SC_REPORT_WARNING("SystemC-AMS",str.str().c_str());
}



void sca_l::matrix_stamps()
{
	//nadd = i
	nadd = add_equation();

	if (!sca_ac_analysis::sca_ac_is_running() && dc_init)
	{
		if (psi0.get() == sca_util::SCA_UNDEFINED)
		{
			//short cut
			B_wr(nadd, p) = 1.0;
			B_wr(nadd, n) = -1.0;

			B_wr(p, nadd) = 1.0;
			B_wr(n, nadd) = -1.0;
		}
		else
		{
			B_wr(nadd, nadd) = -1.0;

			double ph_l = psi0.get() / value.get();

			//test for NaN and infinity
			if (ph_l == std::numeric_limits<double>::infinity() || ph_l
					== -std::numeric_limits<double>::infinity() || ph_l != ph_l)
			{
				// (phi>0)/0 -> infinity
				if (ph_l == ph_l)
				{
					std::ostringstream str;
					str
							<< "The value of the inductance cannot be zero if psi0 is greater zero "
							<< "(results in an infinite initial current) for: "
							<< name();
					SC_REPORT_ERROR("SystemC-AMS",str.str().c_str());
				}

				q(nadd).set_value(0.0);

				B_wr(p,nadd) = 1.0;
				B_wr(n,nadd) = -1.0;


			}
			else
			{
				q(nadd).set_value(ph_l);

				B_wr(p,nadd) = 1.0;
				B_wr(n,nadd) = -1.0;
			}

		}
		add_method(POST_SOLVE, SCA_VMPTR(sca_l::post_solve));
	}
	else
	{
		A(nadd, nadd) = -value.get(); // L*dphi = v

		B_wr(nadd, p) = 1.0;
		B_wr(nadd, n) = -1.0;

		B_wr(p, nadd) = 1.0;
		B_wr(n, nadd) = -1.0;
	}
}


void sca_l::post_solve()
{
	//reset dc stamps
	B_wr(nadd,nadd) = 0.0;
	q(nadd).set_value(0.0);

	B_wr(p,nadd) = 0.0;
	B_wr(n,nadd) = 0.0;


	remove_method(POST_SOLVE, SCA_VMPTR(sca_l::post_solve));

	//set transient stamps
	A(nadd,nadd) = -value.get();  // L*dphi = v

	B_wr(nadd,p) =  1.0;
	B_wr(nadd,n) = -1.0;


	B_wr(p,nadd) =  1.0;
	B_wr(n,nadd) = -1.0;

	request_reinit(1);

	dc_init=false;
}


bool sca_l::trace_init(sca_util::sca_implementation::sca_trace_object_data& data)
{
    data.type=through_value_type;
    data.unit=through_value_unit;

    //trace will be activated after every complete cluster calculation
    //by the synchronization layer
    get_sync_domain()->add_solver_trace(data);
    return true;
}

void sca_l::trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)
{
    sca_core::sca_time ctime = sca_eln::sca_module::get_time();

    double through_value = x(nadd);

	if((id<0)&&(this->trd!=NULL)) //interactive trace
	{
	    this->trd->store_value(through_value);
		return;
	}

    buffer.store_time_stamp(id,ctime,through_value);
}

sca_util::sca_complex sca_l::calculate_ac_result(sca_util::sca_complex* res_vec)
{
	return res_vec[nadd];
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


} //namespace sca_eln

