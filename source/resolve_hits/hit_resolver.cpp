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

#include "hit_resolver.h"

#include <boost/core/ignore_unused.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/sub_range.hpp>

#include "common/algorithm/contains.h"
#include "common/algorithm/sort_uniq_copy.h"
#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/range/adaptor/equal_grouped.h"
#include "common/boost_addenda/range/front.h"
#include "common/c++14/cbegin_cend.h"
#include "resolve_hits/hit_list.h"
#include "resolve_hits/masked_bests_cacher.h"
#include "resolve_hits/scored_hit_arch.h"

#include <map>

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using boost::ignore_unused;
using boost::none;
using boost::range::sort;
using boost::sub_range;
using std::make_pair;
using std::numeric_limits;

// TODO:
//  5. Create scored_arch_proxy (with score and vector of hit indices)
//  6. Change both stores to store that
//  7. Still need to sort out get_arrows_before_starts_of_doms_right_interspersed_with_all_of ?
//  8. Run on all of Dave's hard files. Any problems? Slowest? Any with nested discontigs? Any with interspersed discontigs?
//  9. Consider un-inlining get_best_scored_arch_with_one_of_hits()
// 10. Sort out parsing
// 11. Why is sort allocating?
// 12. Replace sort with in-place insertion during parsing?
// 13. Parse and insert in-place with two, producer/consumer threads?
// 14. Skip inserts of already strictly worse hits?
// 15. Remove strictly worse elements?

/// \brief TODOCUMENT
auto hit_resolver::get_hit_stops_differ_fn(const hit_list &arg_hit_list ///< TODOCUMENT
                                           ) {
	return [&] (const hitidx_t &x, const hitidx_t &y) {
		return ( arg_hit_list[ x ].get_stop_arrow() != arg_hit_list[ y ].get_stop_arrow() );
	};
}


/// \brief TODOCUMENT
///
/// \todo Consider dropping this check it doesn't fire when
///       run under debug mode over some large data set
inline void sanity_check(const scored_arch_proxy &arg_scored_arch_proxy, ///< TODOCUMENT
                         const hit_list          &arg_hits,              ///< TODOCUMENT
                         const hit_vec           &arg_masks              ///< TODOCUMENT
                         ) {
	ignore_unused( arg_scored_arch_proxy, arg_hits, arg_masks );
#ifndef NDEBUG
	for (const auto &arch_hit_idx : arg_scored_arch_proxy) {
		if ( hit_overlaps_with_any_of_hits( arg_hits[ arch_hit_idx ], arg_masks ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("ERROR: Cannot assume that best architecture so far does not clash with masks"));
		}
	}
#endif
}

/// \brief TODOCUMENT
///
/// Find the best possible score and corresponding architecture that fits in the specified region
/// and doesn't overlap with any segments in the specified discontiguous hits
///
/// \todo Could be more efficient at rejecting hits that stop in forbidden regions
scored_arch_proxy hit_resolver::get_best_score_and_arch_of_specified_regions(const hit_vec           &arg_masks,          ///< TODOCUMENT
                                                                             const res_arrow         &arg_start_arrow,    ///< TODOCUMENT. Guaranteed to be at the boundary of a segment in arg_masks, or at the very start if arg_masks is empty
                                                                             const res_arrow         &arg_stop_arrow,     ///< TODOCUMENT. Guaranteed to be at the boundary of a segment in arg_masks, or at the very end   if arg_masks is empty
                                                                             const scored_arch_proxy &arg_best_upto_start ///< TODOCUMENT
                                                                             ) {
	// std::cerr << "Entering get_best_score_and_arch_of_specified_regions() with " << arg_masks.size() << " masks" << std::endl;
	best_scan_arches bests{ arg_stop_arrow.res_before() + 1 };
	if ( arg_start_arrow > start_arrow() && arg_best_upto_start.get_score() > INIT_SCORE ) {
		bests.extend_up_to_arrow( arg_start_arrow - 1 );
		bests.add_best_up_to_arrow( arg_start_arrow, arg_best_upto_start );
	}

	masked_bests_cacher the_masked_bests_cacher = make_masked_bests_cacher(
		the_masked_bests_cache,
		arg_masks,
		hits.get(),
		arg_start_arrow
	);

	const auto indices_of_hits = indices_of_hits_that_stop_in_range( hits, arg_start_arrow, arg_stop_arrow );
	for (const auto &indices_of_hits_with_same_stop : indices_of_hits | equal_grouped( get_hit_stops_differ_fn( hits.get() ) ) ) {
	// for (const auto &hits_with_same_stop : hits_stopping_in_range | equal_grouped( get_hit_stops_differ_fn() ) ) {
		const auto current_arrow = hits.get()[ front( indices_of_hits_with_same_stop ) ].get_stop_arrow();
		sanity_check( bests.get_best_scored_arch_so_far(), hits.get(), arg_masks );
		the_masked_bests_cacher.advance_to_pos_with_best_so_far(
			current_arrow,
			bests.get_best_scored_arch_so_far()
		);

		const auto best_prev_score         = bests.extend_up_to_arrow( current_arrow - 1 );
		const auto best_new_score_and_arch = get_best_scored_arch_with_one_of_hits(
			indices_of_hits_with_same_stop,
			arg_masks,
			arg_start_arrow,
			bests,
			best_prev_score
		);

		// TODOCMENT
		if ( best_new_score_and_arch && best_new_score_and_arch->get_score() > best_prev_score ) {
			bests.add_best_up_to_arrow( current_arrow, *best_new_score_and_arch );
		}
		else {
			// or just... TODOCMENT
			// (this may extend over forbidden regions but that's OK because nothing should try to use those results)
			bests.extend_up_to_arrow( current_arrow );
		}
	}

	sanity_check( bests.get_best_scored_arch_so_far(), hits.get(), arg_masks );
	the_masked_bests_cacher.advance_to_end_with_best_so_far( bests.get_best_scored_arch_so_far() );

	return bests.get_best_scored_arch_so_far();
}

// /// \brief TODOCUMENT
// scored_hit_arch hit_resolver::resolve_step() {
// 	best_scan_arches bests{ max_stop };
// 	const auto score_of_hit   = get_hit_score_fn     ( bests );
// 	const auto hit_score_less = get_hit_score_less_fn( bests );

// 	for (const auto &same_stop_group : hits.get() | equal_grouped( get_hit_stops_differ_fn() ) ) {
// 		const auto stop_arrow      = front( same_stop_group ).get_stop_arrow();
// 		const auto best_prev_score = bests.extend_up_to_arrow( stop_arrow - 1 );

// 		const auto best_new_result_itr = max_element(
// 			same_stop_group,
// 			hit_score_less
// 		);
// 		const auto best_new_score = score_of_hit( *best_new_result_itr );
// 		if ( best_new_score > best_prev_score) {
// 			add_domain_with_its_best( bests, *best_new_result_itr );
// 		}
// 	}
// 	return best_score_and_arch_so_far( bests );
// }


/// \brief TODOCUMENT
hit_resolver::hit_resolver(const hit_list &arg_hits ///< TODOCUMENT
                           ) : hits     ( arg_hits                 ),
                               max_stop ( get_max_stop( arg_hits ) ) {

	constexpr hitidx_t max = numeric_limits<hitidx_t>::max();
	if ( arg_hits.size() + 2 > max ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Number of hits is "
			+ ::std::to_string( arg_hits.size() )
			+ "which is too high, given that the type being used for indexing can only hold a maximum value of "
			+ ::std::to_string( max )
			+ ". You could consider changing the hitidx_t type alias and recompiling."
		));
	}
//	BOOST_LOG_TRIVIAL( warning ) << "get_max_stop is : " << max_stop;
}

/// \brief TODOCUMENT
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

/// \brief TODOCUMENT
scored_hit_arch cath::rslv::resolve_hits(const hit_list &arg_hits ///< TODOCUMENT
                                         ) {
	return detail::hit_resolver{ arg_hits }.resolve();
}


