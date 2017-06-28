/// \file
/// \brief The name_set class definitions

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

#include "name_set.hpp"

#include <iostream>

using std::string;
using std::ostream;

/// \brief Generate a string describing the specified name_set
///
/// \relates name_set
string cath::file::to_string(const name_set &arg_name_set ///< The name_set to describe
                             ) {
	return
		  R"(name_set[name_from_acq:")"
		+ arg_name_set.get_name_from_acq()
		+ R"(")"
		+ (
			arg_name_set.get_specified_id()
			? ( R"(, spec_id:")" + *arg_name_set.get_specified_id() + R"(")" )
			: ""
		)
		+ (
			arg_name_set.get_domain_name_from_regions()
			? ( R"(, dom_id:")" + *arg_name_set.get_domain_name_from_regions() + R"(")" )
			: ""
		)
		+ "]";
}

/// \brief Insert a description of the specified name_set into the specified ostream
///
/// \relates name_set
ostream & cath::file::operator<<(ostream        &arg_os,      ///< The ostream into which the description should be inserted
                                 const name_set &arg_name_set ///< The name_set to describe
                                 ) {
	arg_os << to_string( arg_name_set );
	return arg_os;
}
