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

#include "scored_hit_arch.h"

#include "common/algorithm/transform_build.h"
#include "resolve_hits/hit_list.h"

using namespace cath::common;
using namespace cath::rslv;

/// \brief Make a scored_hit_arch from a scored_arch_proxy and hit_list
///
/// \relates scored_hit_arch
scored_hit_arch cath::rslv::make_scored_hit_arch(const scored_arch_proxy &arg_scored_arch_proxy, ///< The scored_arch_proxy which the scored_hit_arch should copy
                                                 const hit_list          &arg_hit_list           ///< The hit_list to which the scored_arch_proxy refers
                                                 ) {
	return {
		arg_scored_arch_proxy.get_score(),
		hit_arch{
			transform_build<hit_vec>(
				arg_scored_arch_proxy,
				[&] (const hitidx_t &x) {
					return arg_hit_list[ x ];
				}
			)
		}
	};
}

