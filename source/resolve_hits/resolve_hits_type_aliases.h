/// \file
/// \brief The resolve_hits type_aliases header

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

#ifndef RESOLVE_HITS_TYPE_ALIASES_H_INCLUDED
#define RESOLVE_HITS_TYPE_ALIASES_H_INCLUDED

#include <boost/optional/optional_fwd.hpp>

namespace cath { namespace rslv { class hit; } }
namespace cath { namespace rslv { class hit_arch; } }
namespace cath { namespace rslv { class hit_seg; } }
namespace cath { namespace rslv { class res_arrow; } }
namespace cath { namespace rslv { class scored_arch_proxy; } }
namespace cath { namespace rslv { class scored_hit_arch; } }

#include <vector>

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		using hitidx_t = unsigned int;

		/// \brief TODOCUMENT
		using hitidx_hitidx_pair = std::pair<hitidx_t, hitidx_t>;

		/// \brief TODOCUMENT
		using hitidx_vec = std::vector<hitidx_t>;

		/// \brief TODOCUMENT
		using residx_t = unsigned int;

		/// \brief TODOCUMENT
		using resarw_t = residx_t;

		/// \brief TODOCUMENT
		using resscr_t = float;

		/// \brief TODOCUMENT
		constexpr resscr_t INIT_SCORE = 0.0;

		/// \brief TODOCUMENT
		using res_arrow_vec = std::vector<res_arrow>;

		/// \brief TODOCUMENT
		using res_arrow_res_arrow_pair = std::pair<res_arrow, res_arrow>;

		/// \brief TODOCUMENT
		using res_arrow_res_arrow_pair_vec = std::vector<res_arrow_res_arrow_pair>;

		/// \brief TODOCUMENT
		using residx_residx_pair = std::pair<residx_t, residx_t>;

		/// \brief TODOCUMENT;
		using residx_residx_pair_vec = std::vector<residx_residx_pair>;

		/// \brief TODOCUMENT
		using resscr_vec = std::vector<resscr_t>;

		/// \brief TODOCUMENT
		using hit_arch_vec = std::vector<hit_arch>;

		using scored_arch_proxy_opt = boost::optional<scored_arch_proxy>;

		/// \brief TODOCUMENT
		using scored_arch_proxy_vec = std::vector<scored_arch_proxy>;

		/// \brief TODOCUMENT
		using scored_hit_arch_opt = boost::optional<scored_hit_arch>;

		/// \brief TODOCUMENT
		using scored_hit_arch_vec = std::vector<scored_hit_arch>;

		/// \brief TODOCUMENT
		using hit_seg_vec = std::vector<hit_seg>;

		/// \brief TODOCUMENT
		using hit_vec = std::vector<hit>;
	}
}

#endif
