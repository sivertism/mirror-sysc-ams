/*****************************************************************************

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

  sca_traceable_object.cpp - description

  Original Author: Karsten Einwich COSEDA Technologies GmbH Dresden

  Created on: Feb 23, 2016

   SVN Version       :  $Revision$
   SVN last checkin  :  $Date$
   SVN checkin by    :  $Author$
   SVN Id            :  $Id$

 *****************************************************************************/

#include <systemc-ams>
#include "scams/utility/tracing/sca_traceable_object.h"

namespace sca_util
{

void sca_traceable_object::set_unit(const std::string& unit_)
{
	unit=unit_;
}


const std::string& sca_traceable_object::get_unit() const
{
	return unit;
}

void sca_traceable_object::set_unit_prefix(const std::string& prefix)
{
	unit_prefix=prefix;
}


const std::string& sca_traceable_object::get_unit_prefix() const
{
	return unit_prefix;
}

void sca_traceable_object::set_domain(const std::string& domain_)
{
	domain=domain_;
}


const std::string& sca_traceable_object::get_domain() const
{
	return domain;
}


} //namespace sca_util
