/// \file
/// \brief The hit_resolver class definitions

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

#include "hit_resolver.hpp"

#include <boost/core/ignore_unused.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/sub_range.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/sort_uniq_copy.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/adaptor/equal_grouped.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/resolve_hits/algo/masked_bests_cacher.hpp"
#include "cath/resolve_hits/calc_hit_list.hpp"
#include "cath/resolve_hits/resolve/naive_greedy_hit_resolver.hpp"
#include "cath/resolve_hits/scored_hit_arch.hpp"

#include <map>

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;
using namespace cath::seq;

using boost::ignore_unused;
using std::numeric_limits;

// POSSIBLY TODO:
//  * Replace sort with in-place insertion during parsing?
//  * Parse and insert in-place with two, producer/consumer threads?
//  * Skip inserts of already strictly worse hits?
//  * Remove strictly worse elements?
//  * Add options to allow categories of matches (eg input file tying match_ids to categories
//    and options to specify weighting to assign to categories; file could possibly permit regexps for match ID)

/// \brief Return a lambda function that returns whether the hits associated with two
///        specified indices have differing stop points
///
/// This is used to group hits that stop at the same boundary
auto hit_resolver::get_hit_stops_differ_fn(const calc_hit_list &prm_calc_hit_list ///< The list of hits to which the indices will refer
                                           ) {
	return [&] (const hitidx_t &x, const hitidx_t &y) {
		return ( get_stop_arrow( prm_calc_hit_list[ x ] ) != get_stop_arrow( prm_calc_hit_list[ y ] ) );
	};
}


/// \brief Sanity check the best result seen so far doesn't conflict with any of the mask
///
/// \todo Consider dropping this check it doesn't fire when
///       run under debug mode over some large data set
static inline void sanity_check(const scored_arch_proxy &prm_scored_arch_proxy, ///< The architecture (scored_arch_proxy) to check
                                const calc_hit_list     &prm_hits,              ///< The list of hits to which the scored_arch_proxy corresponds
                                const calc_hit_vec      &prm_mask               ///< The mask with which to check for conflicts
                                ) {
	ignore_unused( prm_scored_arch_proxy, prm_hits, prm_mask );
#ifndef NDEBUG
	for (const auto &arch_hit_idx : prm_scored_arch_proxy) {
		if ( hit_overlaps_with_any_of_hits( prm_hits[ arch_hit_idx ], prm_mask ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("ERROR: Cannot assume that best architecture so far does not clash with masks"));
		}
	}
#endif
}

/// \brief Get the best architecture (and score) of the specified regions given
///        the best architecture up to the start point
///
/// This is the main dynamic-programming function.
///
/// Find the best possible score and corresponding architecture that fits in the specified region
/// and doesn't overlap with any segments in the specified discontiguous hits
///
/// \todo Could be more efficient at rejecting hits that stop in forbidden regions
scored_arch_proxy hit_resolver::get_best_score_and_arch_of_specified_regions(const calc_hit_vec      &prm_mask,           ///< The active mask defining any no-go areas.
                                                                             const seq_arrow         &prm_start_arrow,    ///< The point at which to start the dynamic-programming scan. Guaranteed to be at the boundary of a segment in prm_masks, or at the very start if prm_masks is empty
                                                                             const seq_arrow         &prm_stop_arrow,     ///< The point at which to stop the dynamic-programming scan. Guaranteed to be at the boundary of a segment in prm_masks, or at the very end   if prm_masks is empty
                                                                             const scored_arch_proxy &prm_best_upto_start ///< The known-best architecture up to the start point
                                                                             ) {
	// Create and initialise a new best_scan_arches for this particular scan
	// (with size of one more than the stop (just to ensure enough space)
	//  and pre-extended to start with prm_best_upto_start at prm_start_arrow)
	best_scan_arches bests{ prm_stop_arrow.res_before() + 1 };
	if ( prm_start_arrow > start_arrow() && prm_best_upto_start.get_score() > INIT_SCORE ) {
		bests.extend_up_to_arrow( prm_start_arrow - 1 );
		bests.add_best_up_to_arrow( prm_start_arrow, prm_best_upto_start );
	}

	// Create a the_masked_bests_cache that's configured to store results in the_masked_bests_cache
	// whenever it passes a point that may be needed later on
	masked_bests_cacher the_masked_bests_cacher = make_masked_bests_cacher(
		the_masked_bests_cache,
		prm_mask,
		the_dhibs,
		prm_start_arrow
	);

	// Get a range of the indices of the sub-range of the hits that stop within ( prm_start_arrow, prm_stop_arrow ]
	const auto indices_of_hits = indices_of_hits_that_stop_in_range( hits, prm_start_arrow, prm_stop_arrow );

	// Loop over the groups of hits' indices that correspond to hits with the same stop point
	for (const auto &indices_of_hits_with_same_stop : indices_of_hits | equal_grouped( get_hit_stops_differ_fn( hits.get() ) ) ) {
		// Grab the stop point of the hits in this group
		const auto current_arrow = get_stop_arrow( hits.get()[ front( indices_of_hits_with_same_stop ) ] );

		sanity_check( bests.get_best_scored_arch_so_far(), hits.get(), prm_mask );

		// Advance the_masked_bests_cacher to this point (so it stores any results that might be needed later on)
		the_masked_bests_cacher.advance_to_pos_with_best_so_far(
			current_arrow,
			bests.get_best_scored_arch_so_far()
		);

		// In case there was a gap between the previous arrow we considered (or maybe `start_arrow()` / `prm_start_arrow - 1`)
		// and this current arrow, update bests to record that there were no improvements on the previous result
		// up to and including one position before this current_arrow
		// (which will do nothing if there's no gap). Also grab that best previous score.
		//
		// (this may extend over forbidden regions but that's OK because nothing should try to use those results)
		const auto best_prev_score         = bests.extend_up_to_arrow( current_arrow - 1 );

		// See whether it's possible to achieve a better score than that best previous by
		// using one of these hits that stop at current_arrow
		const auto best_new_score_and_arch = get_best_scored_arch_with_one_of_hits(
			indices_of_hits_with_same_stop,
			prm_mask,
			prm_start_arrow,
			bests,
			best_prev_score
		);

		// If we were able to find a result and it has a better score then update bests with that new result
		if ( best_new_score_and_arch && best_new_score_and_arch->get_score() > best_prev_score ) {
			bests.add_best_up_to_arrow( current_arrow, *best_new_score_and_arch );
		}
		// Otherwise the previous result is still the best, so just extend bests to current_arrow
		else {
			bests.extend_up_to_arrow( current_arrow );
		}
	}

	// Sanity check
	sanity_check( bests.get_best_scored_arch_so_far(), hits.get(), prm_mask );

	// Extend the_masked_bests_cacher to the end in case we still need to cache
	// results for use later on
	the_masked_bests_cacher.advance_to_end_with_best_so_far( bests.get_best_scored_arch_so_far() );

	// Return the best result found
	return bests.get_best_scored_arch_so_far();
}

/// \brief Ctor for hit_resolver
hit_resolver::hit_resolver(const calc_hit_list &prm_hits ///< The hits to resolve
                           ) : hits     ( prm_hits                               ),
                               max_stop ( get_max_stop( prm_hits ).value_or( 0 ) ),
                               the_dhibs( prm_hits                               ) {
	constexpr hitidx_t max = numeric_limits<hitidx_t>::max();
	if ( prm_hits.size() + 2 > max ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Number of hits is "
			+ ::std::to_string( prm_hits.size() )
			+ "which is too high, given that the type being used for indexing can only hold a maximum value of "
			+ ::std::to_string( max )
			+ ". You could consider changing the hitidx_t type alias and recompiling."
		));
	}
}

/// \brief Method for resolving hits
///
/// \pre resolve() has never been called before
scored_hit_arch hit_resolver::resolve() {
	return make_scored_hit_arch(
		get_best_score_and_arch_of_specified_regions(
			{},
			start_arrow(),
			arrow_after_res( max_stop ),
			scored_arch_proxy{}
		),
		hits.get()
	);
}

/// \brief The front-end for resolving hits
scored_hit_arch cath::rslv::resolve_hits(const calc_hit_list &prm_hits,        ///< The hits to resolve
                                         const bool          &prm_naive_greedy ///< Whether to use a naive, greedy approach to resolving
                                         ) {
	return prm_naive_greedy ? naive_greedy_resolve_hits( prm_hits )
	                        : detail::hit_resolver{ prm_hits }.resolve();
}


