/// \file
/// \brief The masked_bests_cacher class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHER_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHER_HPP

#include <boost/range/sub_range.hpp>

#include "cath/resolve_hits/algo/masked_bests_cache.hpp"

// clang-format off
namespace cath::rslv { class calc_hit_list; }
namespace cath::rslv::detail { class discont_hits_index_by_start; }
// clang-format on

namespace cath::rslv::detail {

	/// \brief Manage storing best architectures for unmasked-region signatures when the scanning passes
	///        pre-calculated boundaries (at which the results will be needed later on)
	class masked_bests_cacher final {
	private:
		/// \brief A const_iterator type alias for a const_iterator over a res_arrow_vec
		using const_iterator = seq::res_arrow_vec::const_iterator;

		/// \brief (A reference_wrapper to) the masked_bests_cached in which best results should be stored as appropriate
		std::reference_wrapper<masked_bests_cache> cache_ref;

		/// \brief (A reference_wrapper to) the vector of hits that form the mask
		std::reference_wrapper<const calc_hit_vec> masks_ref;

		/// \brief The points at which best results should be stored
		const seq::res_arrow_vec arrows_to_store;

		/// \brief The current position
		const_iterator current_itr;

		void advance_to_itr_with_best_so_far(const const_iterator &,
		                                     const scored_arch_proxy &);

	public:
		masked_bests_cacher(masked_bests_cache &,
		                    const calc_hit_vec &,
		                    seq::res_arrow_vec);

		void advance_to_pos_with_best_so_far(const seq::seq_arrow &,
		                                     const scored_arch_proxy &);

		void advance_to_end_with_best_so_far(const scored_arch_proxy &);
	};

	/// \brief Implementation method to advance to the specified iterator location with the specified best architecture seen so far
	inline void masked_bests_cacher::advance_to_itr_with_best_so_far(const const_iterator    &prm_itr,        ///< The iterator location to which to advance
	                                                                 const scored_arch_proxy &prm_best_so_far ///< The best architecture (scored_arch_proxy) seen so far
	                                                                 ) {
		for (const seq::seq_arrow &the_arrow : boost::sub_range<const seq::res_arrow_vec>( current_itr, prm_itr ) ) {
			store_best_for_masks_up_to_arrow(
				cache_ref.get(),
				prm_best_so_far,
				masks_ref.get(),
				the_arrow
			);
		}
		current_itr = prm_itr;
	}

	/// \brief Ctor from the cache, masks and cache-points
	///
	/// \pre `is_sorted( prm_arrows )`
	inline masked_bests_cacher::masked_bests_cacher(masked_bests_cache  &prm_masked_bests_cache, ///< The cache to which the cacher should store
	                                                const calc_hit_vec  &prm_masks,              ///< The currently-active masks that will define the unmasked-region signatures
	                                                seq::res_arrow_vec   prm_arrows              ///< The points at which to store results in the cache
	                                                ) : cache_ref       ( prm_masked_bests_cache            ),
	                                                    masks_ref       ( prm_masks                         ),
	                                                    arrows_to_store ( std::move( prm_arrows           ) ),
	                                                    current_itr     ( ::std::cbegin( arrows_to_store ) ) {
	}

	/// \brief Advance to the specified position with the specified best architecture (scored_arch_proxy)
	///        performing any cache-stores as appropriate
	inline void masked_bests_cacher::advance_to_pos_with_best_so_far(const seq::seq_arrow    &prm_new_position, ///< The position to which to advance
	                                                                 const scored_arch_proxy &prm_best_so_far   ///< The best architecture (scored_arch_proxy) seen thus far
	                                                                 ) {
		advance_to_itr_with_best_so_far(
			std::find_if(
				current_itr,
				::std::cend( arrows_to_store ),
				[&] (const seq::seq_arrow &x) { return x >= prm_new_position; }
			),
			prm_best_so_far
		);
	}

	/// \brief Advance to the end of the region with the specified best architecture (scored_arch_proxy)
	///        performing any cache-stores as appropriate
	inline void masked_bests_cacher::advance_to_end_with_best_so_far(const scored_arch_proxy &prm_best_so_far ///< The best architecture (scored_arch_proxy) seen thus far
	                                                                 ) {
		advance_to_itr_with_best_so_far(
			::std::cend( arrows_to_store ),
			prm_best_so_far
		);
	}

	masked_bests_cacher make_masked_bests_cacher(masked_bests_cache &,
	                                             const calc_hit_vec &,
	                                             const discont_hits_index_by_start &,
	                                             const seq::seq_arrow &);

	seq::res_arrow_vec get_arrows_before_starts_of_doms_right_interspersed_with_all_of(const calc_hit_vec &,
	                                                                                   const discont_hits_index_by_start &,
	                                                                                   const seq::seq_arrow &);

} // namespace cath::rslv::detail

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_MASKED_BESTS_CACHER_HPP
