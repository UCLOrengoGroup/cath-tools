/// \file
/// \brief The naive_greedy_hit_resolver class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_RESOLVE_NAIVE_GREEDY_HIT_RESOLVER_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_RESOLVE_NAIVE_GREEDY_HIT_RESOLVER_HPP

#include "resolve_hits/calc_hit_list.hpp"

#include <functional>

namespace cath { namespace rslv { class scored_hit_arch; } }

namespace cath {
	namespace rslv {

		scored_hit_arch naive_greedy_resolve_hits(const calc_hit_list &);

	} // namespace rslv
} // namespace cath

#endif
