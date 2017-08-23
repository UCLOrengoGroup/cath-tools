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

#include "hierarchy_group.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <string>

using boost::adaptors::transformed;
using boost::algorithm::join;
using std::ostream;
using std::string;

/// \brief Generate a string describing the specified hierarchy_group
///
/// \relates hierarchy_layer
string cath::clust::to_string(const hierarchy_group &arg_hierarchy_group ///< The hierarchy_group to describe
                              ) {
	return "hierarchy_group["
		+ join(
			arg_hierarchy_group
				| transformed( [] (const hierarchy_value &x) {
					return to_string( x );
				} ),
			", "
		)
		+ "]";
}

/// \brief Insert a description of the specified hierarchy_group into the specified ostream
///
/// \relates hierarchy_layer
ostream & cath::clust::operator<<(ostream               &arg_os,             ///< The ostream into which the description should be inserted
                                  const hierarchy_group &arg_hierarchy_group ///< The hierarchy_group to describe
                                  ) {
	arg_os << to_string( arg_hierarchy_group );
	return arg_os;
}
