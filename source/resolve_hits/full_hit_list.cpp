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

#include <boost/optional.hpp>

#include "common/boost_addenda/range/max_proj_element.hpp"
#include "resolve_hits/full_hit.hpp"

#include <string>

using namespace cath::common;
using namespace cath::rslv;

using boost::make_optional;
using boost::none;
using std::ostream;
using std::string;

/// \brief Generate a string describing the specified full_hit_list
///
/// \relates full_hit_list
string cath::rslv::to_string(const full_hit_list &arg_full_hit_list ///< The full_hit_list to describe
                             ) {
	return "full_hit_list["
		+ ::std::to_string( arg_full_hit_list.size() )
		+ "full_hits]";
}

/// \brief Insert a description of the specified full_hit_list into the specified ostream
///
/// \relates full_hit_list
ostream & cath::rslv::operator<<(ostream             &arg_ostream, ///< The ostream into which the description should be inserted
                                 const full_hit_list &arg_full_hit_list ///< The full_hit_list to describe
                                 ) {
	arg_ostream << to_string( arg_full_hit_list );
	return arg_ostream;
}

/// \brief Get the maximum stop residue of all the full_hits in the specified full_hit_list
///        or none if the full_hit_list is empty
///
/// \relates full_hit_list
residx_opt cath::rslv::get_max_stop(const full_hit_list &arg_full_hit_list ///< The full_hit_list to query
                                    ) {
	// Can't just use make_optional(bool, T) because the second part shouldn't be evaluated if arg_full_hit_list's empty
	return arg_full_hit_list.empty()
		? none
		: make_optional( max_proj(
			arg_full_hit_list,
			std::less<>{},
			&cath::rslv::get_stop_res_arrow
		).res_before() );
}
