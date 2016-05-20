/// \file
/// \brief The masked_bests_cacher class definitions

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

#include "masked_bests_cacher.h"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/upper_bound.hpp>

#include "common/algorithm/copy_build.h"
#include "common/algorithm/sort_uniq_copy.h"
#include "resolve_hits/hit_list.h"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::all_of;
using boost::range::max_element;
using boost::range::upper_bound;
using boost::sub_range;

/// \brief TODOCUMENT
masked_bests_cacher cath::rslv::detail::make_masked_bests_cacher(masked_bests_cache &arg_masked_bests_cache, ///< TODOCUMENT
                                                                 const hit_vec      &arg_masks,              ///< TODOCUMENT
                                                                 const hit_list     &arg_hits,               ///< TODOCUMENT
                                                                 const res_arrow    &arg_start_arrow         ///< TODOCUMENT
                                                                 ) {
	return {
		arg_masked_bests_cache,
		arg_masks,
		get_arrows_before_starts_of_doms_right_interspersed_with_all_of(
			arg_masks,
			arg_hits,
			arg_start_arrow
		)
	};
}

//
// cath|current|3gbnL02/124-195-i5_12,0.042 489.20028051289 3344-3374,3480-3491
// cath|current|1mfaH01/251-347-i5_10,0.013 1057.85904164469 3424-3471,3540-3580

/// \brief TODOCUMENT
///
/// This can be made more efficient by just searching through a separate list of
/// just the discontiguous hits.
///
/// \todo Remove #includes of all_of, filtered and max_element if not used here
// get_arrows_before_starts_of_discontigs_that_fit_with_all_of_and_right_intersperse_with_any_of
res_arrow_vec cath::rslv::detail::get_arrows_before_starts_of_doms_right_interspersed_with_all_of(const hit_vec   &arg_masks,      ///< TODOCUMENT
                                                                                                  const hit_list  &arg_hits,       ///< TODOCUMENT
                                                                                                  const res_arrow &arg_start_arrow ///< TODOCUMENT
                                                                                                  ) {
	// If no masks, just return an empty vector
	if ( arg_masks.empty() ) {
		return {};
	}

	// Otherwise, get the max last stop arrow of any of arg_masks
	const res_arrow &max_stop_arrow = max_element(
		arg_masks,
		hit::get_hit_stop_less()
	)->get_stop_arrow();

	// Find the start of the hits that stop beyond max_stop
	const auto upper_bound_itr = upper_bound(
		arg_hits,
		max_stop_arrow,
		[] (const res_arrow &a, const hit &h) {
			return a < h.get_stop_arrow();
		}
	);

	// const bool is_the_one = ( get_start_res_index( arg_masks.front() ) == 3344 && get_stop_res_index( arg_masks.front() ) == 3491 );
	// if ( is_the_one ) {
	// 	std::cerr << "searching for store arrows for :";
	// 	for (const auto &x : arg_masks ) {
	// 		std::cerr << " " << to_string( x );
	// 	}
	// 	std::cerr << "\n";
	// }

	// Prepare a predicate function that checks whether a hit right intersperses all of arg_masks
	const auto hit_right_intersperses_all_arg_masks = [&] (const hit &x) {
		// if ( is_the_one && x.is_discontig() ) {
			
		// 	const bool answer = (
		// 		x.is_discontig()
		// 		&&
		// 		all_of( arg_masks, [&] (const hit &y) { return second_right_intersperses_first( y, x ); } )
		// 	);
		// 	std::cerr << "\tConsider " << x << " : " << std::boolalpha << answer << "\n";
		// }
		return (
			x.is_discontig()
			&&
			all_of( arg_masks, [&] (const hit &y) { return second_right_intersperses_first( y, x ); } )
		);
	};

	// Return a sorted-uniqued vector of the res_arrows at the start of the hits
	return sort_uniq_copy( copy_build<res_arrow_vec>(
		sub_range<const hit_list>( upper_bound_itr, common::cend( arg_hits ) )
			| filtered   ( hit_right_intersperses_all_arg_masks                     )
			| transformed( [ ] (const hit       &x) { return x.get_start_arrow(); } )
			| filtered   ( [&] (const res_arrow &x) { return x>= arg_start_arrow; } )
	) );
}