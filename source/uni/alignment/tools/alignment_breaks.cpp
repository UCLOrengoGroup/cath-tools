/// \file
/// \brief The alignment_io class definitions

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


#include "alignment_breaks.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/irange.hpp>

#include "alignment/alignment.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/set_union_build.hpp"
#include "common/algorithm/sets_are_disjoint.hpp"
#include "common/boost_addenda/range/adaptor/adjacented.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/size_t_literal.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;

using boost::adaptors::filtered;
using boost::irange;

/// \brief TODOCUMENT
size_vec cath::align::get_alignment_breaks(const alignment &prm_alignment ///< TODOCUMENT
                                           ) {
	// For each alignment index *after* zero, check whether there's any overlap
	// between the entries that are present at the previous index and the entries that are
	// present at this index. If not, add this index to the list.
	//
	// Note: Would prefer to use Boost Range's adjacent_filtered adaptor here but
	// it doesn't do what's required here: it always lets the first pair through,
	// even if they fail the predicate.
	return copy_build<size_vec>(
		irange( 1_z, prm_alignment.length() )
			| filtered (
				[&] (const size_t &curr_idx) {
					// Return whether the set of present entries for the previous index
					// is disjoint with the set of present entries for this index
					return sets_are_disjoint(
						entries_present_at_index( prm_alignment, curr_idx - 1 ),
						entries_present_at_index( prm_alignment, curr_idx     )
					);
				}
			)
	);
}

/// \brief TODOCUMENT
break_pair_validity_and_future cath::align::detail::check_pair(const alignment &prm_alignment, ///< TODOCUMENT
                                                               const size_t    &prm_break_one, ///< TODOCUMENT
                                                               const size_t    &prm_break_two  ///< TODOCUMENT
                                                               ) {
	// LM LR MR
	//  0  0  0  ok; can extend either way
	//  0  0  1  no; can extend right
	/// 0  1  0  no; cannot extend
	/// 0  1  1  no; cannot extend
	//  1  0  0  no; can extend left
	//  1  0  1  no; cannot extend
	/// 1  1  0  no; cannot extend
	/// 1  1  1  no; cannot extend

	const size_vec lefts  = entries_present_at_index( prm_alignment, prm_break_one - 1 );
	const size_vec rights = entries_present_at_index( prm_alignment, prm_break_two     );
	if ( ! sets_are_disjoint( lefts, rights ) ) {
		return { break_pair_validity::BAD, break_pair_future::NEVER_AGAIN };
	}

	const bool left_and_right_are_all = ( lefts.size() + rights.size() == prm_alignment.num_entries() );
	const bool middle_is_empty        = ( prm_break_one == prm_break_two );
	if ( middle_is_empty ) {
		const auto the_future = ( left_and_right_are_all ? break_pair_future::NEVER_AGAIN : break_pair_future::MAYBE_LATER );
		return { break_pair_validity::GOOD, the_future };
	}
	if ( left_and_right_are_all ) {
		return { break_pair_validity::BAD, break_pair_future::NEVER_AGAIN };
	}

	const size_vec middles = entries_present_in_index_range( prm_alignment, prm_break_one, prm_break_two );
	assert( ! middles.empty() );
	if ( ! sets_are_disjoint( lefts, middles ) ) {
		return { break_pair_validity::BAD, break_pair_future::NEVER_AGAIN };
	}
	const auto lefts_and_rights         = set_union_build<size_vec>( lefts,   rights           );
	const auto lefts_middles_and_rights = set_union_build<size_vec>( middles, lefts_and_rights );
	const bool lmr_are_all              = lefts_middles_and_rights.size() == prm_alignment.num_entries();
	const bool middle_conflicts_right = ! sets_are_disjoint( middles, rights );
	return {
		( middle_conflicts_right ? break_pair_validity::BAD       : break_pair_validity::GOOD      ),
		( lmr_are_all            ? break_pair_future::NEVER_AGAIN : break_pair_future::MAYBE_LATER )
	};
}

/// \brief TODOCUMENT
size_size_pair_vec cath::align::get_alignment_break_pairs(const alignment &prm_alignment ///< TODOCUMENT
                                                          ) {
	const auto alignment_breaks     = get_alignment_breaks( prm_alignment );
	const auto num_alignment_breaks = alignment_breaks.size();

	size_size_pair_vec break_pairs;
	break_pairs.reserve( num_alignment_breaks );

	for (const size_t &idx_one : indices( num_alignment_breaks ) ) {
		for (const size_t &idx_two : irange( idx_one, num_alignment_breaks ) ) {
			const size_t &break_one = alignment_breaks[ idx_one ];
			const size_t &break_two = alignment_breaks[ idx_two ];

			const auto pair_result = check_pair( prm_alignment, break_one, break_two );
			if ( pair_result.first  == break_pair_validity::GOOD ) {
				break_pairs.emplace_back( break_one, break_two );
			}
			if ( pair_result.second == break_pair_future::NEVER_AGAIN ) {
				break;
			}
		}
	}
	return break_pairs;
}
