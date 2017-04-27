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

#include "masked_bests_cacher.hpp"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/upper_bound.hpp>

#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_copy.hpp"
#include "exception/out_of_range_exception.hpp"
#include "resolve_hits/algo/discont_hits_index_by_start.hpp"
#include "resolve_hits/calc_hit_list.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;
using namespace cath::seq;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::all_of;
using boost::algorithm::any_of;
using boost::range::max_element;
using boost::range::upper_bound;
using boost::sub_range;

/// \brief Build a masked_bests_cacher to handle the caching of best results in the specified masked_bests_cache
///        whilst scanning from the specified start position with the specified hits and mask
///
/// \relates masked_bests_cacher
masked_bests_cacher cath::rslv::detail::make_masked_bests_cacher(masked_bests_cache                &arg_masked_bests_cache, ///< The cache to which the masked_bests_cacher should store best architectures
                                                                 const calc_hit_vec                &arg_masks,              ///< The currently-active masks that will define the unmasked-region signatures
                                                                 const discont_hits_index_by_start &arg_dhibs,              ///< A discont_hits_index_by_start that allows us to quickly find the discontiguous domains that start within a specified region
                                                                 const seq_arrow                   &arg_start_arrow         ///< The start arrow of the region that will be scanned
                                                                 ) {
	return {
		arg_masked_bests_cache,
		arg_masks,
		get_arrows_before_starts_of_doms_right_interspersed_with_all_of(
			arg_masks,
			arg_dhibs,
			arg_start_arrow
		)
	};
}

/// \brief Calculate the start boundaries of any discontiguous hits that are right interspersed with
///        all of the hits in the mask
///
/// \relates masked_bests_cacher
res_arrow_vec cath::rslv::detail::get_arrows_before_starts_of_doms_right_interspersed_with_all_of(const calc_hit_vec                &arg_masks,      ///< The mask domains, all of which should be right-interspersed by the domains we're looking for
                                                                                                  const discont_hits_index_by_start &arg_dhibs,      ///< A discont_hits_index_by_start that allows us to quickly find the discontiguous domains that start within a specified region
                                                                                                  const seq_arrow                   &arg_start_arrow ///< The start arrow that must precede all of the returned arrows
                                                                                                  ) {
	// If no masks, just return an empty vector
	if ( arg_masks.empty() ) {
		return {};
	}

#ifndef NDEBUG
	// If in debug mode, check all hits in the mask are discontiguous, or throw otherwise
	if ( ! all_of( arg_masks, [] (const calc_hit &x) { return is_discontig( x ); } ) ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Mask should only contain discontiguous domains"));
	}
#endif

	// Grab the latest start point of the mask's hits' interiors
	const seq_arrow &max_first_seg_stop = get_stop_of_first_segment( *max_element(
		arg_masks,
		calc_hit::get_hit_first_seg_stop_less()
	));

	// Grab the earliest stop point of the mask's hits' interiors
	const seq_arrow &min_last_seg_start = get_start_of_last_segment( *max_element(
		arg_masks,
		calc_hit::get_hit_last_seg_start_less()
	));

#ifndef NDEBUG
	// If in debug mode, check that the latest gap start is strictly before the earliest gap stop, or throw otherwise
	if ( max_first_seg_stop >= min_last_seg_start ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Mask's hits' interiors should have a non-empty common region"));
	}
#endif

	// Use the discont_hits_index_by_start to find the discontiguous domains that start
	// in this region - ask it to return the relevant indices in its own index
	const auto index_indices_of_disconts_in_range = arg_dhibs.get_index_indices_of_disconts_in_range(
		max_first_seg_stop,
		min_last_seg_start
	);

	// Prepare a predicate function that checks whether a hit suitable intersperses all of arg_masks
	//
	// This can't just be `all_of(... second_right_intersperses_first(...) ...)` because it must
	// handle situations like:
	// ~~~
	// c <--->                           <--->
	// a      <--->           <--->
	//                  ^
	//                  |
	//                  v
	// b                 <--->     <--->
	// ~~~~
	// in which b doesn't right intersperse a but the arrow point should still be recognised
	// as a relevant start in the mask c+a
	const auto hit_suitable_intersperses_arg_masks = [&] (const calc_hit &x) {
		return (
			is_discontig( x )
			&&
			all_of( arg_masks, [&] (const calc_hit &y) { return second_right_or_inside_intersperses_first( y, x ); } )
			&&
			any_of( arg_masks, [&] (const calc_hit &y) { return second_right_intersperses_first( y, x ); } )
		);
	};

	// Return a sorted-uniqued vector of the res_arrows at the start of the hits
	return sort_uniq_copy( copy_build<res_arrow_vec>(
		index_indices_of_disconts_in_range
			| transformed( [&] (const size_t    &x) { return arg_dhibs.get_discont_hit_of_index_index( x ); } )
			| filtered   ( hit_suitable_intersperses_arg_masks                      )
			| transformed( [ ] (const calc_hit  &x) { return get_start_arrow( x ); } )
			| filtered   ( [&] (const seq_arrow &x) { return x>= arg_start_arrow; } )
	) );
}
