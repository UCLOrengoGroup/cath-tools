/// \file
/// \brief The data_option definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "data_option.h"

#include "exception/invalid_argument_exception.h"

using namespace cath::common;
using namespace cath::opts::detail;
using namespace std;

/// \brief Simple insertion operator for data_option
///
/// \relates data_option
ostream & cath::opts::detail::operator<<(ostream                &arg_os,    ///< The ostream to which to output the data_option
                                         const data_option &arg_data_option ///< The data_option to output
                                         ) {
	switch (arg_data_option) {
		case ( data_option::PATH   ) : { arg_os << "data_option::PATH"   ; break ; }
		case ( data_option::PREFIX ) : { arg_os << "data_option::PREFIX" ; break ; }
		case ( data_option::SUFFIX ) : { arg_os << "data_option::SUFFIX" ; break ; }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of data_option not recognised whilst inserting into an ostream"));
			break; // Superfluous, post-throw break statement to appease Eclipse's syntax highlighter
		}
	}
	return arg_os;
}
