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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHE_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHE_HPP

#include <unordered_map>

#include "cath/common/hash/hash_value_combine.hpp"
#include "cath/common/type_traits.hpp"
#include "cath/resolve_hits/algo/scored_arch_proxy.hpp"
#include "cath/resolve_hits/calc_hit.hpp"
#include "cath/seq/seq_seg.hpp"

namespace cath::rslv {

	namespace detail {

		/// \brief Provide function operator that hashes seq_seg_vec so that they can be used as
		///        keys in an unordered map
		struct seq_seg_vec_hasher final {

			/// \brief Hash function for seq_seg_vec
			size_t operator()(const seq::seq_seg_vec &prm_seq_seg_vec ///< The seq_seg_vec to hash
			                  ) const {
				size_t seed = 0;
				for (const seq::seq_seg &the_hit : prm_seq_seg_vec) {
					for (const auto &index : { the_hit.get_start_arrow().get_index(),
					                           the_hit.get_stop_arrow ().get_index(), } ) {
						using index_type = common::remove_cvref_t< decltype( index ) >;
						common::hash_value_combine( seed, std::hash< index_type >{}( index ) );
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
		inline seq::seq_seg_vec get_unmasked_regions_before_arrow(const calc_hit_vec   &prm_hits, ///< The hits defining the mask. These must be non-overlapping but may be unsorted.
		                                                          const seq::seq_arrow &prm_arrow ///< The point at which to stop
		                                                          ) {
			// Get a sorted copy of prm_hits's segments
			const auto seq_segs = get_start_sorted_seq_segs( prm_hits );

			// Prepare the working data: a seq_seg_vec to populate and a seq_arrow at the end of the most-recently-handled seq_seg
			seq::seq_seg_vec results;
			auto prev_stop = seq::start_arrow();

			// Loop over the mask segments
			for (const seq::seq_seg &the_seq_seg : seq_segs) {

				// If this mask segment starts after the stop arrow, then break out of the loop
				if ( the_seq_seg.get_start_arrow() >= prm_arrow ) {
					break;
				}
				// Else if this mask segment starts *strictly* after the previous stop, add a record for the gap
				if ( the_seq_seg.get_start_arrow() >  prev_stop ) {
					results.emplace_back( prev_stop, the_seq_seg.get_start_arrow() );
				}
				// Update the prev_stop to this segment's stop
				prev_stop = the_seq_seg.get_stop_arrow();
			}

			// If the stop point is *strictly* after the previously handled segment's stop, add a record for the gap
			// (this happens in all cases except those where there is a mask segment stopping-at or straddling prm_arrow)
			if ( prm_arrow > prev_stop ) {
				results.emplace_back( prev_stop, prm_arrow );
			}

			return results;
		}

	} // namespace detail

	/// \brief Store the best scored_arch_proxy for a given unmasked pattern
	class masked_bests_cache final {
	private:
		/// \brief The unordered map (ie hash-map to store the optimal architecture (scored_arch_proxy) for a given set of unmasked regions (seq_seg_vec))
		std::unordered_map<seq::seq_seg_vec, scored_arch_proxy, detail::seq_seg_vec_hasher> store;

	public:
		const scored_arch_proxy & get_best_for_unmasked(const seq::seq_seg_vec &) const;
		void store_best_for_unmasked(seq::seq_seg_vec &&,
		                             const scored_arch_proxy &);
	};

	/// \brief Get the optimum architecture (scored_arch_proxy) for the specified signature of unmasked regions
	inline const scored_arch_proxy & masked_bests_cache::get_best_for_unmasked(const seq::seq_seg_vec &prm_unmasked ///< The set of unmasked regions for which the optimum architecture is required
	                                                                           ) const {
		return store.at( prm_unmasked );
	}

	/// \brief Store the optimum architecture for a signature of unmasked regions
	inline void masked_bests_cache::store_best_for_unmasked(seq::seq_seg_vec        &&prm_unmasked,              ///< The set of unmasked regions for which the optimum architecture is to be stored
	                                                        const scored_arch_proxy  &prm_best_scored_arch_proxy ///< The optimum architecture (scored_arch_proxy) to store
	                                                        ) {
		store.emplace(
			std::move( prm_unmasked ),
			prm_best_scored_arch_proxy
		);
	}

	/// \brief Get the optimum architecture (scored_arch_proxy) from the specified masked_bests_cache
	///        for the signature of regions unmasked by the specified mask up to the specified point
	///
	/// \relates masked_bests_cache
	inline const scored_arch_proxy & get_best_for_masks_up_to_arrow(const masked_bests_cache &prm_masked_bests_cache, ///< The masked_bests_cache to query
	                                                                const calc_hit_vec       &prm_mask_hits,          ///< The mask that defines the unmasked regions for which the architecture is optimal
	                                                                const seq::seq_arrow     &prm_stop_arrow          ///< The stop boundary at which the signature of unmasked regions should stop
	                                                                ) {
		return prm_masked_bests_cache.get_best_for_unmasked(
			detail::get_unmasked_regions_before_arrow( prm_mask_hits, prm_stop_arrow )
		);
	}

	/// \brief Store the optimum architecture (scored_arch_proxy) in the specified masked_bests_cache
	///        for the signature of regions unmasked by the specified mask up to the specified point
	///
	/// \relates masked_bests_cache
	///
	///     Check       : nothing in arch - discontigs should overrun current_arrow
	///     Check       : nothing previously stored there
	inline void store_best_for_masks_up_to_arrow(masked_bests_cache      &prm_masked_bests_cache,     ///< The masked_bests_cache in which to store the optimum architecture
	                                             const scored_arch_proxy &prm_best_scored_arch_proxy, ///< The optimum architecture for the unmasked regions implied by the other arguments
	                                             const calc_hit_vec      &prm_mask_hits,              ///< The mask that defines the unmasked regions for which the architecture is optimal
	                                             const seq::seq_arrow    &prm_stop_arrow              ///< The stop boundary at which the signature of unmasked regions should stop
	                                             ) {
		prm_masked_bests_cache.store_best_for_unmasked(
			detail::get_unmasked_regions_before_arrow( prm_mask_hits, prm_stop_arrow ),
			prm_best_scored_arch_proxy
		);
	}

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHE_HPP
