/// \file
/// \brief The get_sorting_scores class definitions

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

#include "get_sorting_scores.hpp"

#include "common/algorithm/sorted_indices.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/invert_permutation.hpp"

using namespace cath;
using namespace cath::common;

/// \brief Get the ordering rank of each entry
///
///        (ie the position in the resulting size_vec that corresponds to the "first"
///         element has value 0)
size_vec cath::clust::get_sorting_scores(const id_of_str_bidirnl &arg_name_ider, ///< The name/ID lookup for the entries (used in ascending sorting)
                                         const doub_vec          &arg_props      ///< The properties on which to sort the entries (used in ascending sorting)
                                         ) {
	if ( arg_props.size() != arg_name_ider.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get sorting_scores for a name_ider and properties list that have different sizes"));
	}
	return common::invert_permutation( common::sorted_indices(
		arg_name_ider.size(),
		[&] (const size_t &x, const size_t &y) {
			return (
				tie( arg_props[ x ], arg_name_ider.get_name_of_id( x ) )
				<
				tie( arg_props[ y ], arg_name_ider.get_name_of_id( y ) )
			);
		}
	) );
}

/// \brief Get the ordering rank of each entry
///
///        (ie the position in the resulting size_vec that corresponds to the "first"
///         element has value 0)
size_vec cath::clust::get_sorting_scores(const id_of_str_bidirnl &arg_name_ider ///< The name/ID lookup for the entries (used in ascending sorting)
                                         ) {
	return common::invert_permutation( common::sorted_indices(
		arg_name_ider.size(),
		[&] (const size_t &x, const size_t &y) {
			return (
				tie( arg_name_ider.get_name_of_id( x ) )
				<
				tie( arg_name_ider.get_name_of_id( y ) )
			);
		}
	) );
}
