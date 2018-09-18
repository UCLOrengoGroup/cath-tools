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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_RESOLVE_HIT_RESOLVER_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_RESOLVE_HIT_RESOLVER_HPP

#include <boost/optional.hpp>
#include <boost/range/sub_range.hpp>

#include "resolve_hits/algo/best_scan_arches.hpp"
#include "resolve_hits/algo/discont_hits_index_by_start.hpp"
#include "resolve_hits/algo/masked_bests_cache.hpp"
#include "resolve_hits/calc_hit_list.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

#include <functional>

namespace cath { namespace rslv { class scored_hit_arch; } }

namespace cath {
	namespace rslv {

		namespace detail {

			/// \brief Store the data used in the resolve hits algorithm and provide that algorithm
			///
			/// Use this via calls to cath::rslv::resolve_hits()
			///
			/// \todo Separate out these two different responsibilities
			///
			/// \todo What happens if client calls resolve() twice?
			///
			/// Like all cath-resolve-hits code, this assumes simple residue numbering
			/// and is hence unsuitable for use with raw PDB residue numbers.
			class hit_resolver final {
			private:
				static auto get_hit_stops_differ_fn(const calc_hit_list &);

				static void update_best_if_hit_improves(scored_arch_proxy_opt &,
				                                        const resscr_t &,
			                                            const hitidx_t &,
				                                        const scored_arch_proxy &,
				                                        const resscr_t &);

				scored_arch_proxy_opt get_best_scored_arch_with_one_of_hits(const boost::sub_range<boost::integer_range<hitidx_t>> &,
				                                                            const calc_hit_vec &,
				                                                            const seq::seq_arrow &,
				                                                            const best_scan_arches &,
				                                                            const resscr_t &);

				/// \brief A reference to the hits to be resolved
				std::reference_wrapper<const calc_hit_list> hits;

				/// \brief The maximum stop of any of the hits
				seq::residx_t max_stop;

				/// \brief A cache of the best architectures seen for specific unmasked-region signatures
				masked_bests_cache the_masked_bests_cache;

				/// \brief An index of the discontiguous domains by their start residues
				///        which enables much faster searches for right-interspersers of masks
				///
				/// This is initialised on construction and isn't modified after that
				discont_hits_index_by_start the_dhibs;

				scored_arch_proxy get_best_score_and_arch_of_specified_regions(const calc_hit_vec &,
				                                                               const seq::seq_arrow &,
				                                                               const seq::seq_arrow &,
				                                                               const scored_arch_proxy &);

			public:
				explicit hit_resolver(const calc_hit_list &);

				scored_hit_arch resolve();
			};

			/// \brief Update the scored_arch_proxy_opt if the specified hit improves
			///        on the best result seen so far
			///
			/// To update, the new score must beat:
			///  * prm_best_so_far->get_score() if prm_best_so_far is set
			///  * prm_score_to_beat            otherwise
			inline void hit_resolver::update_best_if_hit_improves(scored_arch_proxy_opt   &prm_best_so_far,         ///< The best architecture and score seen so far
			                                                      const resscr_t          &prm_hit_score,           ///< The score of the new hit
			                                                      const hitidx_t          &prm_hit_index,           ///< The index of the new hit
			                                                      const scored_arch_proxy &prm_best_hit_complement, ///< The best architecture to complement this new hit
			                                                      const resscr_t          &prm_score_to_beat        ///< The base score that any result must beat
			                                                      ) {
				const resscr_t this_score = prm_hit_score + prm_best_hit_complement.get_score();
				const bool improves = prm_best_so_far ? ( this_score > prm_best_so_far->get_score() )
				                                      : ( this_score > prm_score_to_beat            );
				if ( improves ) {
					prm_best_so_far = add_hit_copy( prm_best_hit_complement, prm_hit_score, prm_hit_index );
				}
			}

			/// \brief Get the best architecture that can be achieved using one of the specified hits
			///        that all stop at the same boundary
			///
			/// \todo Consider getting prm_start_arrow from prm_bests
			///
			/// Note: Not using max_element because that'd probably call get_complex_hit_score() twice for many elements
			inline scored_arch_proxy_opt hit_resolver::get_best_scored_arch_with_one_of_hits(const boost::sub_range<boost::integer_range<hitidx_t>> &prm_hit_indices,  ///< The indices of the hits to consider
			                                                                                 const calc_hit_vec                                     &prm_mask,         ///< The active mask defining no-go regions
			                                                                                 const seq::seq_arrow                                   &prm_start_arrow,  ///< The start point of the current scan (from which prm_bests should have results (up to the one place before the stop of these hits))
			                                                                                                                                                           ///< Guaranteed to be at the boundary of a segment in prm_masks, or at the very start if prm_masks is empty
			                                                                                 const best_scan_arches                                 &prm_bests,        ///< The history of best-seen architectures so far in this layer of dynamic programming. This is setup to handle the current forbidden prm_masks.
			                                                                                 const resscr_t                                         &prm_score_to_beat ///< The existing score to beat. This is setup to handle the current forbidden prm_masks.
			                                                                                 ) {
				// Loop over each of the hits that stop at prm_current_arrow
				scored_arch_proxy_opt best_so_far;
				for (const auto &hit_index : prm_hit_indices) {
					const auto &the_hit = hits.get()[ hit_index ];

					// If this hit clashes with the forbidden regions marked out by prm_mask,
					// then it can't be used so skip to the next one
					if ( hit_overlaps_with_any_of_hits( the_hit, prm_mask ) ) {
						continue;
					}

					// If the hit is contiguous, then the score/arch is just the sum of:
					//  * the score/arch of this hit itself
					//  * the best score/arch before the start
					//
					// (and we know that the local best_scan_arches will include an entry for the hit's start_arrow because
					//  the hit must be within this scan because the hit is contiguous and prm_start_arrow is never within
					//  the middle of a free stretch)
					//
					// If that's an improvement on best_so_far/prm_score_to_beat then update best_so_far
					if ( ! is_discontig( the_hit ) ) {
						update_best_if_hit_improves(
							best_so_far,
							the_hit.get_score(),
							hit_index,
							prm_bests.get_best_scored_arch_up_to_arrow( get_start_arrow( the_hit ) ),
							prm_score_to_beat
						);
					}
					// Else if the hit is discontiguous...
					else {
						const auto &hit_start = get_start_arrow( the_hit );

						const scored_arch_proxy &best_hit_complement =
							// If the discontiguous hit is within this region
							//   (ie hit_start >= prm_start_arrow;
							//    we already know it must stop before the end of scan, else we wouldn't be considering it)
							//
							// Then the score/arch is the sum of:
							//  * the score/arch of this hit itself
							//  * the best score/arch that can be found by recursing into a new layer of dynamic-programming:
							//    [ up to the end of the hit that doesn't clash with prm_mask or the_hit;
							//      starting from the start of this hit, by using the (already calculated) best score/arch up to that point ]
							( hit_start >= prm_start_arrow ) ? get_best_score_and_arch_of_specified_regions(
								prm_mask + the_hit,
								get_stop_of_first_segment( the_hit ),
								get_start_of_last_segment( the_hit ),
								prm_bests.get_best_scored_arch_up_to_arrow( hit_start )
							) :

							// Else this is a discontiguous hit that starts before prm_start_arrow
							//
							// Then the score/arch is the sum of:
							//  * the score/arch of this hit itself
							//  * the best score/arch that can be found by recursing into a new layer of dynamic programming:
							//    [ up to the end of the hit that doesn't clash with prm_mask or the_hit;
							//      starting from start_arrow again, by using the previously recorded partial pickup of
							//      best score/arch up to there that doesn't clash with (prm_mask + the_hit) ]
							get_best_score_and_arch_of_specified_regions(
								prm_mask + the_hit,
								prm_start_arrow,
								get_start_of_last_segment( the_hit ),
								get_best_for_masks_up_to_arrow( the_masked_bests_cache, prm_mask + the_hit, prm_start_arrow )
							);

						// If that's an improvement on best_so_far/prm_score_to_beat then update best_so_far
						update_best_if_hit_improves(
							best_so_far,
							the_hit.get_score(),
							hit_index,
							best_hit_complement,
							prm_score_to_beat
						);
					}
				}

				// Return result
				return best_so_far;
			}


		} // namespace detail

		scored_hit_arch resolve_hits(const calc_hit_list &,
		                             const bool &);

	} // namespace rslv
} // namespace cath

#endif
