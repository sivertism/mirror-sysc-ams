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

  sca_eln_isource.cpp - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden / COSEDA Technologies

  Created on: 12.11.2009

   SVN Version       :  $Revision: 1944 $
   SVN last checkin  :  $Date: 2016-03-11 14:28:43 +0100 (Fri, 11 Mar 2016) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_eln_isource.cpp 1944 2016-03-11 13:28:43Z karsten $

 *****************************************************************************/

/*****************************************************************************/


#include "systemc-ams"
#include "scams/predefined_moc/eln/sca_eln_isource.h"
#include "scams/impl/util/data_types/sca_method_object.h"
#include "scams/impl/util/data_types/sca_method_list.h"
#include "scams/impl/core/sca_solver_base.h"
#include "scams/impl/predefined_moc/conservative/sca_con_interactive_trace_data.h"

#include "scams/impl/analysis/ac/sca_ac_domain_globals.h"
#include "scams/impl/analysis/ac/sca_ac_domain_db.h"
#include "scams/impl/analysis/ac/sca_ac_domain_eq.h"

using namespace sca_ac_analysis;
using namespace sca_implementation;


namespace sca_eln
{

sca_isource::sca_isource(sc_core::sc_module_name, double init_value_ ,
		double offset_ , double amplitude_ , double frequency_ , double phase_,
		const sca_core::sca_time& delay_ , double ac_amplitude_ ,
		double ac_phase_ , double ac_noise_amplitude_ ) :
	p("p"), n("n"), init_value("init_value", init_value_), offset("offset",
			offset_), amplitude("amplitude", amplitude_), frequency(
			"frequency", frequency_), phase("phase",phase_), delay("delay", delay_), ac_amplitude(
			"ac_amplitude", ac_amplitude_),
			ac_phase("ac_phase", ac_phase_), ac_noise_amplitude(
					"ac_noise_amplitude", ac_noise_amplitude_)
{
    through_value_available = true;
    through_value_type      = "I";
    through_value_unit      = "A";

    i_value=0.0;

    nadd = -1;
}

const char* sca_isource::kind() const
{
	return "sca_eln::sca_isource";
}

void sca_isource::trace( sc_core::sc_trace_file* tf ) const
{
	std::ostringstream str;
	str << "sc_trace of sca_isource module not supported for module: " << this->name();
	SC_REPORT_WARNING("SystemC-AMS",str.str().c_str());
}



double sca_isource::i_t()
{
    if (get_time() < delay)
    {
    	i_value=init_value;
    	return i_value;
    }

    i_value = offset + amplitude * sin(2.0 * M_PI * frequency * (sca_get_time()
    			- delay.get().to_seconds()) + phase);
    return i_value;
}


void sca_isource::matrix_stamps()
{
    nadd = add_equation();

    B_wr(p, nadd) =  1.0;
    B_wr(n, nadd) = -1.0;

    if( sca_ac_analysis::sca_ac_is_running() )
    {
    	if(sca_ac_analysis::sca_ac_noise_is_running())
    	{
            B_wr(nadd, nadd) =  1.0;
            A(nadd, nadd) =  0.0;

        	if(get_ac_database().is_initialized() && get_ac_database().set_src_value())
        	{
        		long id;
        		id=get_ac_database().get_start_of_add_equations(this->get_sync_domain());

        		if(id>=0)
        		{
        			get_ac_database().get_equations().get_noise_src_names()(id+nadd)=this->name();
        			get_ac_database().get_equations().get_Bgnoise()(id+nadd)=ac_noise_amplitude.get();
        		}
        	}
    	}
    	else
    	{
    		sca_util::sca_complex phase_cmpl(0, (M_PI/180) * ac_phase.get());

    		phase_cmpl = exp(-phase_cmpl);

    		B_wr(nadd, nadd) =  phase_cmpl.real();
    		A(nadd, nadd) =  phase_cmpl.imag()/sca_ac_analysis::sca_ac_w();

    		q(nadd).set_value(-ac_amplitude.get());
    	}
    }
    else
    {
        B_wr(nadd, nadd) =  1.0;
        A(nadd, nadd) =  0.0;

        q(nadd).sub_element(SCA_MPTR(sca_isource::i_t), this);
    }

}



bool sca_isource::trace_init(sca_util::sca_implementation::sca_trace_object_data& data)
{
    data.type=through_value_type;
    data.unit=through_value_unit;

    //trace will be activated after every complete cluster calculation
    //by the synchronization layer
    get_sync_domain()->add_solver_trace(data);
    return true;
}

void sca_isource::trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)
{
    sca_core::sca_time ctime = sca_eln::sca_module::get_time();

    double through_value = i_value;

	if((id<0)&&(this->trd!=NULL)) //interactive trace
	{
	    this->trd->store_value(through_value);
		return;
	}

    buffer.store_time_stamp(id,ctime,through_value);
}

sca_util::sca_complex sca_isource::calculate_ac_result(sca_util::sca_complex* res_vec)
{
	return std::polar(ac_amplitude.get(),(M_PI/180.0) * ac_phase.get());
}



const double& sca_isource::get_typed_trace_value() const
{
	if(this->trd==NULL)
	{
		this->trd=new sca_core::sca_implementation::sca_con_interactive_trace_data(this);

	}

	return this->trd->get_value();
}


const std::string& sca_isource::get_trace_value() const
{
	std::ostringstream str;
	str << this->get_typed_trace_value();

	static std::string trace_value_string;
	trace_value_string=str.str().c_str();
	return trace_value_string;
}





bool sca_isource::register_trace_callback(sca_trace_callback cb,void* cb_arg)
{
	if(this->trd==NULL)
	{
		this->trd=new sca_core::sca_implementation::sca_con_interactive_trace_data(this);
	}


	this->trd->register_trace_callback(cb,cb_arg);

	return true;

}


} //namespace sca_eln
