/// \file
/// \brief The hit_seg_boundary_fns class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_TRIM_HIT_SEG_BOUNDARY_FNS_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_TRIM_HIT_SEG_BOUNDARY_FNS_H

#include "resolve_hits/resolve_hits_type_aliases.h"

namespace cath { namespace rslv { class trim_spec; } }

namespace cath {
	namespace rslv {

		namespace detail {

			/// \brief Represent which boundary is being requested
			enum class boundary_wanted : bool {
				START, ///< The start boundary is wanted
				STOP   ///< The stop boundary is wanted
			};

			res_arrow_opt get_boundary_impl(const boundary_wanted &,
			                                const hit_seg &,
			                                const hit_seg_vec &,
			                                const trim_spec &);
		} // namespace detail

		res_arrow_opt get_start_boundary(const hit_seg &,
		                                 const hit_seg_vec &,
		                                 const trim_spec &);

		res_arrow_opt get_stop_boundary(const hit_seg &,
		                                const hit_seg_vec &,
		                                const trim_spec &);

		seg_boundary_pair get_boundary_pair(const hit_seg &,
		                                    const hit_seg_vec &,
		                                    const trim_spec &);

		res_arrow_opt calc_resolved_boundary(const hit_seg &,
		                                     const hit_seg &,
		                                     const trim_spec &);

	} // namespace rslv
} // namespace cath

#endif
