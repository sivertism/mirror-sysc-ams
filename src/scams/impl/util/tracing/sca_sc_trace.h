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

  sca_sc_trace.h - description

  Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

  Created on: 09.12.2009

   SVN Version       :  $Revision: 1265 $
   SVN last checkin  :  $Date: 2011-11-19 21:43:31 +0100 (Sat, 19 Nov 2011) $
   SVN checkin by    :  $Author: karsten $
   SVN Id            :  $Id: sca_sc_trace.h 1265 2011-11-19 20:43:31Z karsten $

 *****************************************************************************/

/*****************************************************************************/


#ifndef SCA_SC_TRACE_H_
#define SCA_SC_TRACE_H_

namespace sca_util
{
namespace sca_implementation
{

class sca_sc_trace_base : public sca_util::sca_traceable_object
{
public:
	sca_sc_trace_base();

	bool trace_init(sca_trace_object_data& data);
	virtual void trace(long id,sca_trace_buffer& buffer) = 0;
	virtual const sc_core::sc_event& get_trigger_event() = 0;

	virtual void set_type_info(sca_trace_object_data& data)=0;

	sca_trace_object_data* data;

	virtual  ~sca_sc_trace_base();

	bool* terminated;

};



template<typename TIF>
class sca_sc_type_extractor
{
public:
	sca_sc_type_extractor(sca_trace_object_data& data){}
};

template<typename T>
class sca_sc_type_extractor<sc_core::sc_port<sc_core::sc_signal_in_if<T> > >
{
public:
	sca_sc_type_extractor(sca_trace_object_data& data)
	{
		data.set_type_info<T>();
	}
};

template<typename T>
class sca_sc_type_extractor<sc_core::sc_port<sc_core::sc_signal_inout_if<T> > >
{
public:
	sca_sc_type_extractor(sca_trace_object_data& data)
	{
		data.set_type_info<T>();
	}
};


template<typename T>
class sca_sc_type_extractor<sc_core::sc_signal_in_if<T> >
{
public:
	sca_sc_type_extractor(sca_trace_object_data& data)
	{
		data.set_type_info<T>();
	}
};





template<class TIF>
class sca_sc_trace : public sca_sc_trace_base
{

public:
	sca_sc_trace(const TIF& sig_port) : inp(&sig_port)
	{
	}

private:
	void trace(long id, sca_trace_buffer& buffer)
	{
		buffer.store_time_stamp(id, sc_core::sc_time_stamp(),inp->read());
	}

	const sc_core::sc_event& get_trigger_event()
	{
		return inp->default_event();
	}

	void set_type_info(sca_trace_object_data& data)
	{
		sca_sc_type_extractor<TIF> tp(data);
	}

	const TIF* inp;
};

template<class T>
class sca_sc_trace<sc_core::sc_port<T> > : public sca_sc_trace_base
{

public:
	sca_sc_trace(const sc_core::sc_port<T>& sig_port) : inp(&sig_port)
	{
	}

private:
	void trace(long id, sca_trace_buffer& buffer)
	{
		buffer.store_time_stamp(id, sc_core::sc_time_stamp(),(*inp)->read());
	}

	const sc_core::sc_event& get_trigger_event()
	{
		return (*inp)->default_event();
	}

	void set_type_info(sca_trace_object_data& data)
	{
		sca_sc_type_extractor<sc_core::sc_port<T> > tp(data);
	}

	const sc_core::sc_port<T>* inp;
};






}
}


#endif /* SCA_SC_TRACE_H_ */
