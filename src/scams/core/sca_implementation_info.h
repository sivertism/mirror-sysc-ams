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

 sca_implementation_info.h - functions for accessing implementation info

 Original Author: Karsten Einwich Fraunhofer IIS/EAS Dresden

 Created on: 08.03.2009

 SVN Version       :  $Revision: 1984 $
 SVN last checkin  :  $Date: 2016-04-04 13:10:24 +0200 (Mon, 04 Apr 2016) $
 SVN checkin by    :  $Author: karsten $
 SVN Id            :  $Id: sca_implementation_info.h 1984 2016-04-04 11:10:24Z karsten $

 *****************************************************************************/
/*
 * LRM clause 6.4.
 * Implementation information
 */

/*****************************************************************************/

#ifndef SCA_IMPLEMENTATION_INFO_H_
#define SCA_IMPLEMENTATION_INFO_H_


namespace sca_core
{
#define IEEE_16661_SYSTEMC_AMS 201601L

#define SCA_VERSION_MAJOR        2
#define SCA_VERSION_MINOR        1
#define SCA_VERSION_PATCH        0
#define SCA_VERSION_ORGINATOR    "COSEDA"
#define SCA_VERSION_RELEASE_DATE "20160404"
#define SCA_VERSION_PRERELEASE   "RELEASE"
#define SCA_IS_PRERELEASE        0
#define SCA_VERSION              "2.1.0-COSEDA"
#define SCA_COPYRIGHT            "        Copyright (c) 2010-2014  by Fraunhofer-Gesellschaft IIS/EAS\n" \
	                             "        Copyright (c) 2015-2016  by COSEDA Technologies GmbH\n" \
                                 "        Licensed under the Apache License, Version 2.0"



extern const unsigned int sca_version_major;
extern const unsigned int sca_version_minor;
extern const unsigned int sca_version_patch;
extern const std::string  sca_version_orginator;
extern const std::string  sca_version_release_date;
extern const std::string  sca_version_prerelease;
extern const bool         sca_is_prerelease;
extern const std::string  sca_version_string;
extern const std::string  sca_copyright_string;

const char* sca_copyright();
const char* sca_version();
const char* sca_release();

} // namespace sca_core


#endif /* SCA_IMPLEMENTATION_INFO_H_ */
