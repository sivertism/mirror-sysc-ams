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

 sca_eln_tdf_isink.h - electrical linear net converter from current to a tdf signal

 Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

 Created on: 08.03.2009

 SVN Version       :  $Revision: 1523 $
 SVN last checkin  :  $Date: 2013-02-17 21:36:57 +0100 (Sun, 17 Feb 2013) $ (UTC)
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: sca_eln_tdf_isink.h 1523 2013-02-17 20:36:57Z karsten $

 *****************************************************************************/
/*
 * LRM clause 4.3.1.25.
 * The class sca_eln::sca_tdf::sca_isink shall implement a primitive module for
 * the ELN MoC that realizes a scaled conversion from an ELN current to a TDF
 * output signal. The value of the current flowing through the primitive from
 * terminal p to terminal n shall be scaled with coefficient scale and written
 * to a TDF output port outp. The primitive shall contribute the following
 * equation to the equation system:
 *
 *        v(p,n)(t) = 0
 *
 * where v(p,n)(t) is the voltage across terminals p and n.
 */

/*****************************************************************************/

#ifndef SCA_ELN_TDF_ISINK_H_
#define SCA_ELN_TDF_ISINK_H_

namespace sca_eln
{

namespace sca_tdf
{

//    class sca_isink : public implementation-derived-from sca_core::sca_module,
//                      protected sca_util::sca_traceable_object
class sca_isink: public    sca_eln::sca_module,
		         public sca_util::sca_traceable_object
{
public:
	sca_eln::sca_terminal p;
	sca_eln::sca_terminal n;

	::sca_tdf::sca_out<double> outp;

	sca_core::sca_parameter<double> scale;

	virtual const char* kind() const;

	explicit sca_isink(sc_core::sc_module_name, double scale_ = 1.0);


	//begin implementation specific

	private:
		 virtual void matrix_stamps();
		 void assign_result();

		 long nadd;

		 bool trace_init(sca_util::sca_implementation::sca_trace_object_data& data);
		 void trace(long id,sca_util::sca_implementation::sca_trace_buffer& buffer);
		 sca_util::sca_complex calculate_ac_result(sca_util::sca_complex* res_vec);

		 //sc trace of sc_core::sc_object to prevent clang warning due overloaded
		 //sca_traceable_object function
		 void trace( sc_core::sc_trace_file* tf ) const;

	//end implementation specific
};

} // namespace sc_tdf

typedef sca_eln::sca_tdf::sca_isink sca_tdf_isink;

} // namespace sca_eln

#endif /* SCA_ELN_TDF_ISINK_H_ */
