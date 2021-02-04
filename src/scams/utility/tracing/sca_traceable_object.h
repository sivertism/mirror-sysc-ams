/*****************************************************************************

    Copyright 2010
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

 sca_traceable_object.h - base class for all objects which are traceable

 Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

 Created on: 08.03.2009

 SVN Version       :  $Revision: 1914 $
 SVN last checkin  :  $Date: 2016-02-23 19:06:06 +0100 (Tue, 23 Feb 2016) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: sca_traceable_object.h 1914 2016-02-23 18:06:06Z karsten $

 *****************************************************************************/
/*
 * LRM clause 6.1.1.8.
 * The class sca_util::sca_traceable_object shall be the base class for all
 * objects, which can be traced.
 */

/*****************************************************************************/

#ifndef SCA_TRACEABLE_OBJECT_H_
#define SCA_TRACEABLE_OBJECT_H_

namespace sca_util
{
namespace sca_implementation
{
class sca_trace_object_data;
class sca_trace_buffer;
}

class sca_traceable_object
{


public:
	virtual bool trace_init(sca_util::sca_implementation::sca_trace_object_data& data)=0;
	virtual void trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer)=0;


	typedef void (*sca_trace_callback)(void*);

	/** introspection interface */

	/** adds callback for introspection */
	virtual bool register_trace_callback(sca_trace_callback,void*)
	{ return false;}

	/** returns current value */
	virtual const std::string& get_trace_value() const { return empty_string;}

	virtual bool force_value(const std::string&) {return false;}

	virtual void release_value(){}

	//required for Microsoft Visual Compiler 2008
	sca_traceable_object(){}
	virtual ~sca_traceable_object(){}

	virtual sca_util::sca_complex calculate_ac_result(sca_util::sca_complex* res_vec)
	{
	   	  return sca_util::sca_complex(0.0,0.0);
	}


	void set_unit(const std::string& unit);
	const std::string& get_unit() const;

	void set_unit_prefix(const std::string& prefix);
	const std::string& get_unit_prefix() const;

	void set_domain(const std::string& domain);
	const std::string& get_domain() const;

protected:

	/*implementation-defined*/

	const std::string empty_string;

private:


    std::string unit;
    std::string unit_prefix;
    std::string domain;

};

} // namespace sca_util

#endif /* SCA_TRACEABLE_OBJECT_H_ */
