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

#include <vector>

namespace cath { namespace rslv { class hit_seg; } }
namespace cath { namespace rslv { class hit; } }
namespace cath { namespace rslv { class res_arrow; } }
namespace cath { namespace rslv { class scored_arch_proxy; } }

namespace cath {
	namespace rslv {

		/// \brief Type alias for a vector of hit_seg objects
		using hit_seg_vec = std::vector<hit_seg>;

		/// \brief Type alias for the type to be used to index hits
		using hitidx_t = unsigned int;

		/// \brief Type alias for a vector of hitidx_t values
		using hitidx_vec = std::vector<hitidx_t>;

		/// \brief Type alias for a vector of hit objects
		using hit_vec = std::vector<hit>;

		/// \brief Type alias for a vector of res_arrow objects
		using res_arrow_vec = std::vector<res_arrow>;

		/// \brief Type alias for the type to be used to index residues
		using residx_t = unsigned int;

		/// \brief Type alias for a pair of res_arrow and hit index
		using res_arr_idx_pair = std::pair<res_arrow, hitidx_t>;

		/// \brief Type alias for a vector of pairs of res_arrow and hitidx_t
		using res_arr_idx_pair_vec = std::vector<res_arr_idx_pair>;

		/// \brief Type alias for the type to be used for indexing residue boundaries
		///        (within a res_arrow object)
		using resarw_t = residx_t;

		/// \brief Type alias for a pair of residue indices
		using residx_residx_pair = std::pair<residx_t, residx_t>;

		/// \brief Type alias for a vector of pairs of residue indices
		using residx_residx_pair_vec = std::vector<residx_residx_pair>;

		/// \brief Type alias for the type to be used for hits' scores
		using resscr_t = float;

		/// \brief Type alias for an optional scored_arch_proxy object
		using scored_arch_proxy_opt = boost::optional<scored_arch_proxy>;

		/// \brief Type alias for a vector of scored_arch_proxy objects
		using scored_arch_proxy_vec = std::vector<scored_arch_proxy>;

		/// \brief The initial score before any hits have been added
		constexpr resscr_t INIT_SCORE = 0.0;

	}
}

#endif
