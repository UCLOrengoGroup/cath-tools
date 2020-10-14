/// \file
/// \brief The hierarchy class definitions

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

#include "hierarchy_value.hpp"

#include <iostream>
#include <string>

using ::std::ostream;
using ::std::string;

/// \brief Generate a string describing the specified hierarchy_value
///
/// \relates hierarchy_layer
string cath::clust::to_string(const hierarchy_value &prm_hierarchy_value ///< The hierarchy_value to describe
                              ) {
	using ::std::to_string;
	return ( prm_hierarchy_value.get_type() == hierarchy_ref::CLUSTER )
		? "deeper_group_" + to_string( prm_hierarchy_value.get_index() )
		: to_string( prm_hierarchy_value.get_index() );
}

/// \brief Insert a description of the specified hierarchy_value into the specified ostream
///
/// \relates hierarchy_layer
ostream & cath::clust::operator<<(ostream               &prm_os,             ///< The ostream into which the description should be inserted
                                  const hierarchy_value &prm_hierarchy_value ///< The hierarchy_value to describe
                                  ) {
	prm_os << to_string( prm_hierarchy_value );
	return prm_os;
}
