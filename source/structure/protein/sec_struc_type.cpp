/// \file
/// \brief The sec_struc_type definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
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

#include "sec_struc_type.hpp"

#include <boost/throw_exception.hpp>

#include "exception/invalid_argument_exception.hpp"

#include <istream>
#include <string>

using namespace cath;
using namespace cath::common;
using namespace std;

/// \brief Simple extraction operator for sec_struc_type that expects H or S
///
/// \relates sec_struc_type
istream & cath::operator>>(istream        &arg_istream,       ///< The istream from which to extract the sec_struc_type
                           sec_struc_type &arg_sec_struc_type ///< The sec_struc_type to populate
                           ) {
	string input_string;
	arg_istream >> input_string;
	if (input_string == "H") {
		arg_sec_struc_type = sec_struc_type::ALPHA_HELIX;
	}
	else if (input_string == "S") {
		arg_sec_struc_type = sec_struc_type::BETA_STRAND;
	}
	else {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to recognise sec_struc type " + input_string));
	}
	return arg_istream;
}

/// \brief Simple insertion operator for sec_struc_type that outputs H or S
///
/// \relates sec_struc_type
ostream & cath::operator<<(ostream              &arg_ostream,       ///< The ostream to which to output the sec_struc_type
                           const sec_struc_type &arg_sec_struc_type ///< The sec_struc_type to output
                           ) {
	switch(arg_sec_struc_type) {
		case(sec_struc_type::ALPHA_HELIX) : {
			arg_ostream << "H";
			break;
		}
		case(sec_struc_type::BETA_STRAND) : {
			arg_ostream << "S";
			break;
		}
		case(sec_struc_type::COIL) : {
			arg_ostream << " ";
			break;
		}
	}
	return arg_ostream;
}
