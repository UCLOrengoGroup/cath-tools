/// \file
/// \brief The full_hit_list class definitions

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

#include "full_hit_list.hpp"

#include <boost/range/algorithm/equal.hpp>

#include "cath/common/boost_addenda/range/max_proj_element.hpp"
#include "cath/common/optional/make_optional_if.hpp"
#include "cath/resolve_hits/full_hit.hpp"

#include <string>

using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::boost::range::equal;
using ::std::make_optional;
using ::std::nullopt;
using ::std::ostream;
using ::std::string;

/// \brief Generate a string describing the specified full_hit_list
///
/// \relates full_hit_list
string cath::rslv::to_string(const full_hit_list &prm_full_hit_list ///< The full_hit_list to describe
                             ) {
	return "full_hit_list["
		+ ::std::to_string( prm_full_hit_list.size() )
		+ " full_hits]";
}

/// \brief Insert a description of the specified full_hit_list into the specified ostream
///
/// \relates full_hit_list
ostream & cath::rslv::operator<<(ostream             &prm_ostream,      ///< The ostream into which the description should be inserted
                                 const full_hit_list &prm_full_hit_list ///< The full_hit_list to describe
                                 ) {
	prm_ostream << to_string( prm_full_hit_list );
	return prm_ostream;
}

/// \brief Return whether the two specified full_hit_lists are identical
///
/// \relates full_hit_list
bool cath::rslv::operator==(const full_hit_list &prm_lhs, ///< The first  full_hit_list to compare
                            const full_hit_list &prm_rhs  ///< The second full_hit_list to compare
                            ) {
	return equal( prm_lhs, prm_rhs );
}

/// \brief Get the maximum stop residue of all the full_hits in the specified full_hit_list
///        or nullopt if the full_hit_list is empty
///
/// \relates full_hit_list
residx_opt cath::rslv::get_max_stop(const full_hit_list &prm_full_hit_list ///< The full_hit_list to query
                                    ) {
	return if_then_optional(
		! prm_full_hit_list.empty(),
		max_proj(
			prm_full_hit_list,
			std::less<>{},
			&cath::rslv::get_stop_res_arrow
		).res_before()
	);
}
