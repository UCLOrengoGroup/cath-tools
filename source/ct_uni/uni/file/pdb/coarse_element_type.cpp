/// \file
/// \brief The coarse_element_type class definitions

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

#include "coarse_element_type.hpp"

#include "common/exception/invalid_argument_exception.hpp"

using namespace cath::common;

using std::ostream;
using std::string;

/// \brief Generate a string describing the specified coarse_element_type
///
/// \relates coarse_element_type
string cath::file::to_string(const coarse_element_type &prm_element ///< The coarse_element_type to describe
                             ) {
	switch ( prm_element ) {
		case ( coarse_element_type::CARBON       ) : { return "carbon"       ; }
		case ( coarse_element_type::CARBON_ALPHA ) : { return "carbon_alpha" ; }
		case ( coarse_element_type::CARBON_BETA  ) : { return "carbon_beta"  ; }
		case ( coarse_element_type::NITROGEN     ) : { return "nitrogen"     ; }
		case ( coarse_element_type::OXYGEN       ) : { return "oxygen"       ; }
		case ( coarse_element_type::NON_CORE     ) : { return "non_core"     ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of coarse_element_type not recognised whilst converting to_string()"));
}

/// \brief Insert a description of the specified coarse_element_type into the specified ostream
///
/// \relates coarse_element_type
ostream & cath::file::operator<<(ostream                   &prm_os,     ///< The ostream into which the description should be inserted
                                 const coarse_element_type &prm_element ///< The coarse_element_type to describe
                                 ) {
	prm_os << to_string( prm_element );
	return prm_os;
}
