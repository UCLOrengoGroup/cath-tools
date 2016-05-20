/// \file
/// \brief The hit_resolver class header

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

#ifndef HIT_RESOLVER_H_INCLUDED
#define HIT_RESOLVER_H_INCLUDED

#include <boost/optional.hpp>
#include <boost/range/sub_range.hpp>

#include "resolve_hits/best_scan_arches.h"
#include "resolve_hits/hit_list.h"
#include "resolve_hits/masked_bests_cache.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

#include <functional>
//
namespace cath { namespace rslv { class best_scan_arches; } }
namespace cath { namespace rslv { class hit_list; } }

namespace cath {
	namespace rslv {

		namespace detail {

			/// \brief TODOCUMENT
			class hit_resolver final {
			private:
				static auto get_hit_stops_differ_fn(const hit_list &);

				static void update_best_if_hit_improves(scored_arch_proxy_opt &,
				                                        const resscr_t &,
			                                            const hitidx_t &,
				                                        const scored_arch_proxy &,
				                                        const resscr_t &);

				scored_arch_proxy_opt get_best_scored_arch_with_one_of_hits(const boost::sub_range<boost::integer_range<hitidx_t>> &,
				                                                            const hit_vec &,
				                                                            const res_arrow &,
				                                                            const best_scan_arches &,
				                                                            const resscr_t &);

				/// \brief TODOCUMENT
				std::reference_wrapper<const hit_list> hits;

				/// \brief TODOCUMENT
				residx_t max_stop;

				/// \brief TODOCUMENT
				masked_bests_cache the_masked_bests_cache;

				res_arrow_vec get_arrows_before_starts_of_doms_right_interspersed_with_all_of(const hit_vec &);

				scored_arch_proxy get_best_score_and_arch_of_specified_regions(const hit_vec &,
				                                                             const res_arrow &,
				                                                             const res_arrow &,
				                                                             const scored_arch_proxy &);

			public:
				explicit hit_resolver(const hit_list &);

				scored_hit_arch resolve();
			};

			/// \brief TODOCUMENT
			inline void hit_resolver::update_best_if_hit_improves(scored_arch_proxy_opt   &arg_best_so_far,         ///< TODOCUMENT
			                                                      const resscr_t          &arg_hit_score,           ///< TODOCUMENT
			                                                      const hitidx_t          &arg_hit_index,           ///< TODOCUMENT
			                                                      const scored_arch_proxy &arg_best_hit_complement, ///< TODOCUMENT
			                                                      const resscr_t          &arg_score_to_beat        ///< TODOCUMENT
			                                                      ) {
				
				const resscr_t this_score = arg_hit_score + arg_best_hit_complement.get_score();
				const bool improves = arg_best_so_far ? ( this_score > arg_best_so_far->get_score() )
				                                      : ( this_score > arg_score_to_beat            );
				if ( improves ) {
					arg_best_so_far = add_hit_copy( arg_best_hit_complement, arg_hit_score, arg_hit_index );
				}
			}

			/// \brief TODOCUMENT
			///
			/// \todo Consider getting arg_start_arrow from arg_bests
			///
			/// Note: Not using max_element because that'd probably call get_complex_hit_score() twice for many elements
			inline scored_arch_proxy_opt hit_resolver::get_best_scored_arch_with_one_of_hits(const boost::sub_range<boost::integer_range<hitidx_t>> &arg_hit_indices,  ///< TODOCUMENT
			                                                                                 const hit_vec                                          &arg_masks,        ///< TODOCUMENT
			                                                                                 const res_arrow                                        &arg_start_arrow,  ///< The start point of the current scan (from which arg_bests should have results (up to the one place before the stop of these hits))
			                                                                                                                                                           ///< Guaranteed to be at the boundary of a segment in arg_masks, or at the very start if arg_masks is empty
			                                                                                 const best_scan_arches                                 &arg_bests,        ///< TODOCUMENT. This is setup to handle the current forbidden arg_masks.
			                                                                                 const resscr_t                                         &arg_score_to_beat ///< TODOCUMENT. This is setup to handle the current forbidden arg_masks.
			                                                                                 ) {
				// Loop over each of the hits that stop at arg_current_arrow
				scored_arch_proxy_opt best_so_far;
				for (const auto &hit_index : arg_hit_indices) {
					const auto &the_hit = hits.get()[ hit_index ];

					// If this hit clashes with the forbidden regions marked out by arg_masks,
					// then it can't be used so skip
					if ( hit_overlaps_with_any_of_hits( the_hit, arg_masks ) ) {
						continue;
					}

					// If the hit is contiguous, then the score/arch is just the sum of:
					//  * the score/arch of this hit itself
					//  * the best score/arch before the start
					// (The local best_scan_arches must include entry for the hit's start_arrow because
					//  the hit is contiguous and arg_start_arrow is never in never within the middle of a free stretch)
					//
					// If that's an improvement on best_so_far/arg_score_to_beat then update best_so_far
					if ( ! the_hit.is_discontig() ) {
						update_best_if_hit_improves(
							best_so_far,
							the_hit.get_score(),
							hit_index,
							arg_bests.get_best_scored_arch_up_to_arrow( the_hit.get_start_arrow() ),
							arg_score_to_beat
						);
					}
					// Else if the hit is discontiguous...
					else {
						const auto &hit_start = the_hit.get_start_arrow();
						const scored_arch_proxy &best_hit_complement =
							// If ***TODOCUMENT*** (discontiguous, within this region)
							//
							// TODOCUMENT... the sum of:
							//  * the score/arch of this hit itself
							//  * the best score/arch that can be found up to the end of the hit that doesn't clash with arg_masks or the_hit
							//    (starting from the start of this hit, by using the (already calculated) best score/arch up to that point)
							( hit_start >= arg_start_arrow ) ? get_best_score_and_arch_of_specified_regions(
								arg_masks + the_hit,
								get_stop_of_first_segment( the_hit ),
								get_start_of_last_segment( the_hit ),
								arg_bests.get_best_scored_arch_up_to_arrow( hit_start )
							) :

							// Else ***TODOCUMENT***
							//
							// TODOCUMENT... the sum of:
							//  * the score/arch of this hit itself
							//  * the best score/arch that can be found up to the end of the hit that doesn't clash with arg_masks or the_hit
							//    (starting from start_arrow again, by using the previously recorded partial pickup of best score/arch up to there that
							//     doesn't clash with (arg_masks + the_hit)
							get_best_score_and_arch_of_specified_regions(
								arg_masks + the_hit,
								arg_start_arrow,
								get_start_of_last_segment( the_hit ),
								get_best_for_masks_up_to_arrow( the_masked_bests_cache, arg_masks + the_hit, arg_start_arrow )
							);

						// If that's an improvement on best_so_far/arg_score_to_beat then update best_so_far
						update_best_if_hit_improves(
							best_so_far,
							the_hit.get_score(),
							hit_index,
							best_hit_complement,
							arg_score_to_beat
						);
					}
				}

				// Return result
				return best_so_far;
			}


		}

		scored_hit_arch resolve_hits(const hit_list &);


	}
}

#endif
