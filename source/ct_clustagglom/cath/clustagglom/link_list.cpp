/// \file
/// \brief The link_list class definitions

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

#include "link_list.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <string>

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::std::string;

/// \brief Generate a string describing the specified link_list
///
/// \relates links
string cath::clust::link_list_string(const link_list &prm_links, ///< The link_list to describe
                                     const size_t    &prm_base   ///< The base node from which the links emanate
                                     ) {
	using ::std::to_string;
	return to_string( prm_base ) + " -> " + join(
		prm_links
			| transformed( [] (const link &x) {
				return
					  to_string( x.node   )
					+ " ["
					+ to_string( x.dissim )
					+ "]";
			} ),
		", "
	);
}
