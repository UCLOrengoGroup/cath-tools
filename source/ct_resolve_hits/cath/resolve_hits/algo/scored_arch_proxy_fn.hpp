/// \file
/// \brief The scored_arch_proxy functions header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_SCORED_ARCH_PROXY_FN_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_SCORED_ARCH_PROXY_FN_HPP

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/resolve_hits/algo/scored_arch_proxy.hpp"
#include "cath/resolve_hits/calc_hit_list.hpp"

namespace cath::rslv {

	/// \brief Return whether the hit with the specified index overlaps with any of the hits in the specified scored_arch_proxy,
	///        which is tied to the specified calc_hits_list with the specified index 
	///
	/// \relates scored_arch_proxy
	inline bool overlaps_with(const scored_arch_proxy &prm_scored_arch_proxy, ///< The scored_arch_proxy to look for overlaps with
	                          const calc_hit          &prm_hit,               ///< The hit to look for overlaps with
	                          const calc_hit_list     &prm_calc_hit_list      ///< The calc_hit_list to which the scored_arch_proxy is tied
	                          ) {
		return boost::algorithm::any_of(
			prm_scored_arch_proxy
				| boost::adaptors::transformed( [&] (const hitidx_t &x) { return prm_calc_hit_list[ x ]; } ),
			[&] (const calc_hit &x) {
				return are_overlapping( x, prm_hit );
			}
		);
	}

	/// \brief Return whether the specified hit overlaps with any of the hits in the specified scored_arch_proxy,
	///        which is tied to the specified calc_hits_list with the specified index
	///
	/// \relates scored_arch_proxy
	inline bool overlaps_with(const scored_arch_proxy &prm_scored_arch_proxy, ///< The scored_arch_proxy to look for overlaps with
	                          const hitidx_t          &prm_hit_index,         ///< The index of the hit to look for overlaps with
	                          const calc_hit_list     &prm_calc_hit_list      ///< The calc_hit_list to which the scored_arch_proxy is tied
	                          ) {
		return overlaps_with(
			prm_scored_arch_proxy,
			prm_calc_hit_list[ prm_hit_index ],
			prm_calc_hit_list
		);
	}

	/// \brief Check whether the hit with the specified index overlaps with any of the hits in the specified scored_arch_proxy,
	///        (which is tied to the specified calc_hits_list with the specified index)
	///        and if not, add the hit to the proxy
	///
	/// \relates scored_arch_proxy
	inline void add_hit_if_does_not_overlap(scored_arch_proxy   &prm_scored_arch_proxy, ///< The scored_arch_proxy to which the hit should potentially be added
	                                        const resscr_t      &prm_score,             ///< The score associated with the hit to add
	                                        const hitidx_t      &prm_hit_index,         ///< The index of the hit to add
	                                        const calc_hit_list &prm_calc_hit_list      ///< The calc_hit_list to which the scored_arch_proxy is tied
	                                        ) {
		if ( ! overlaps_with( prm_scored_arch_proxy, prm_hit_index, prm_calc_hit_list ) ) {
			prm_scored_arch_proxy.add_hit( prm_score, prm_hit_index );
		}
	}

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_ALGO_SCORED_ARCH_PROXY_FN_HPP
