/// \file
/// \brief The scored_hit_arch class definitions

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

#include "scored_hit_arch.hpp"

#include "common/algorithm/transform_build.hpp"
#include "resolve_hits/calc_hit_list.hpp"

using namespace cath::common;
using namespace cath::rslv;

/// \brief Make a scored_hit_arch from a scored_arch_proxy and calc_hit_list
///
/// \relates scored_hit_arch
scored_hit_arch cath::rslv::make_scored_hit_arch(const scored_arch_proxy &arg_scored_arch_proxy, ///< The scored_arch_proxy which the scored_hit_arch should copy
                                                 const calc_hit_list     &arg_calc_hit_list      ///< The calc_hit_list to which the scored_arch_proxy refers
                                                 ) {
	return {
		arg_scored_arch_proxy.get_score(),
		hit_arch{
			transform_build<calc_hit_vec>(
				arg_scored_arch_proxy,
				[&] (const hitidx_t &x) {
					return arg_calc_hit_list[ x ];
				}
			)
		}
	};
}

/// \brief Get a list of the full hits from the specified full_hit_list that correspond to
///        the hits in the specified scored_hit_arch
///
/// \relates scored_hit_arch
full_hit_list cath::rslv::get_full_hits_of_hit_arch(const scored_hit_arch &arg_scored_hit_arch, ///< The scored_hit_arch whose full hits should be extracted
                                                    const full_hit_list   &arg_full_hits        ///< The full_hit_list associated with the scored_hit_arch, from which the full_hits should be extracted
                                                    ) {
	return get_full_hits_of_hit_arch( arg_scored_hit_arch.get_arch(), arg_full_hits );
}
