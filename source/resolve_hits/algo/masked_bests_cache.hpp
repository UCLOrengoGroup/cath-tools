/// \file
/// \brief The masked_bests_cache class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHE_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHE_H

#include "resolve_hits/algo/scored_arch_proxy.hpp"
#include "resolve_hits/calc_hit.hpp"
#include "resolve_hits/hit_seg.hpp"

#include <unordered_map>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Provide function operator that hashes hit_seg_vec so that they can be used as
			///        keys in an unordered map
			struct hit_seg_vec_hasher final {

				/// \brief Hash function for hit_seg_vec
				size_t operator()(const hit_seg_vec &arg_hit_seg_vec ///< The hit_seg_vec to hash
				                  ) const {
					size_t seed = 0;
					for (const hit_seg &the_hit : arg_hit_seg_vec) {
						for (const auto &index : { the_hit.get_start_arrow().get_index(),
						                           the_hit.get_stop_arrow ().get_index(), } ) {
							using index_type = std::decay_t< decltype( index ) >;
							seed ^= std::hash< index_type >{}( index  ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
						}
					}
					return seed;
				}
			};

			/// \brief Build a list of the regions between zero and the specified arrow
			///        that aren't masked by the specified hits.
			///
			/// \brief The specified hits must be non-overlapping.
			///
			/// Note: this excludes any zero-length regions left by the mask, which means that
			///       it can give identical results for different hit_vecs
			inline hit_seg_vec get_unmasked_regions_before_arrow(const calc_hit_vec &arg_hits, ///< The hits defining the mask. These must be non-overlapping but may be unsorted.
			                                                     const res_arrow    &arg_arrow ///< The point at which to stop
			                                                     ) {
				// Get a sorted copy of arg_hits's segments
				const auto hit_segs = get_start_sorted_hit_segs( arg_hits );

				// Prepare the working data: a hit_seg_vec to populate and a res_arrow at the end of the most-recently-handled hit_seg
				hit_seg_vec results;
				auto prev_stop = start_arrow();

				// Loop over the mask segments
				for (const hit_seg &the_hit_seg : hit_segs) {

					// If this mask segment starts after the stop arrow, then break out of the loop
					if ( the_hit_seg.get_start_arrow() >= arg_arrow ) {
						break;
					}
					// Else if this mask segment starts *strictly* after the previous stop, add a record for the gap
					if ( the_hit_seg.get_start_arrow() >  prev_stop ) {
						results.emplace_back( prev_stop, the_hit_seg.get_start_arrow() );
					}
					// Update the prev_stop to this segment's stop
					prev_stop = the_hit_seg.get_stop_arrow();
				}

				// If the stop point is *strictly* after the previously handled segment's stop, add a record for the gap
				// (this happens in all cases except those where there is a mask segment stopping-at or straddling arg_arrow)
				if ( arg_arrow > prev_stop ) {
					results.emplace_back( prev_stop, arg_arrow );
				}

				return results;
			}
		} // namespace detail

		/// \brief Store the best scored_arch_proxy for a given unmasked pattern
		class masked_bests_cache final {
		private:
			/// \brief The unordered map (ie hash-map to store the optimal architecture (scored_arch_proxy) for a given set of unmasked regions (hit_seg_vec))
			std::unordered_map<hit_seg_vec, scored_arch_proxy, detail::hit_seg_vec_hasher> store;

		public:
			const scored_arch_proxy & get_best_for_unmasked(const hit_seg_vec &) const;
			void store_best_for_unmasked(const hit_seg_vec &&,
			                             const scored_arch_proxy &);
		};

		/// \brief Get the optimum architecture (scored_arch_proxy) for the specified signature of unmasked regions
		inline const scored_arch_proxy & masked_bests_cache::get_best_for_unmasked(const hit_seg_vec &arg_unmasked ///< The set of unmasked regions for which the optimum architecture is required
		                                                                           ) const {
			return store.at( arg_unmasked );
		}

		/// \brief Store the optimum architecture for a signature of unmasked regions
		inline void masked_bests_cache::store_best_for_unmasked(const hit_seg_vec       &&arg_unmasked,              ///< The set of unmasked regions for which the optimum architecture is to be stored
		                                                        const scored_arch_proxy  &arg_best_scored_arch_proxy ///< The optimum architecture (scored_arch_proxy) to store
		                                                        ) {
			store.emplace(
				std::move( arg_unmasked ),
				std::move( arg_best_scored_arch_proxy )
			);
		}

		/// \brief Get the optimum architecture (scored_arch_proxy) from the specified masked_bests_cache
		///        for the signature of regions unmasked by the specified mask up to the specified point
		///
		/// \relates masked_bests_cache
		inline const scored_arch_proxy & get_best_for_masks_up_to_arrow(const masked_bests_cache &arg_masked_bests_cache, ///< The masked_bests_cache to query
		                                                                const calc_hit_vec       &arg_mask_hits,          ///< The mask that defines the unmasked regions for which the architecture is optimal
		                                                                const res_arrow          &arg_stop_arrow          ///< The stop boundary at which the signature of unmasked regions should stop
		                                                                ) {
			return arg_masked_bests_cache.get_best_for_unmasked(
				detail::get_unmasked_regions_before_arrow( arg_mask_hits, arg_stop_arrow )
			);
		}

		/// \brief Store the optimum architecture (scored_arch_proxy) in the specified masked_bests_cache
		///        for the signature of regions unmasked by the specified mask up to the specified point
		///
		/// \relates masked_bests_cache
		///
		///     Check       : nothing in arch - discontigs should overrun current_arrow
		///     Check       : nothing previously stored there
		inline void store_best_for_masks_up_to_arrow(masked_bests_cache      &arg_masked_bests_cache,     ///< The masked_bests_cache in which to store the optimum architecture
		                                             const scored_arch_proxy &arg_best_scored_arch_proxy, ///< The optimum architecture for the unmasked regions implied by the other arguments
		                                             const calc_hit_vec      &arg_mask_hits,              ///< The mask that defines the unmasked regions for which the architecture is optimal
		                                             const res_arrow         &arg_stop_arrow              ///< The stop boundary at which the signature of unmasked regions should stop
		                                             ) {
			arg_masked_bests_cache.store_best_for_unmasked(
				detail::get_unmasked_regions_before_arrow( arg_mask_hits, arg_stop_arrow ),
				arg_best_scored_arch_proxy
			);
		}

	} // namespace rslv
} // namespace cath

#endif
