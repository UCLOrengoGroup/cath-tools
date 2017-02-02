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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_RESOLVE_HITS_TYPE_ALIASES_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_RESOLVE_HITS_TYPE_ALIASES_H

#include <boost/optional/optional_fwd.hpp>

#include <vector>

namespace cath { namespace rslv { class calc_hit; } }
namespace cath { namespace rslv { class crh_segment_spec; } }
namespace cath { namespace rslv { class full_hit; } }
namespace cath { namespace rslv { class hit_seg; } }
namespace cath { namespace rslv { class res_arrow; } }
namespace cath { namespace rslv { class scored_arch_proxy; } }
namespace cath { namespace rslv { struct alnd_rgn; } }
namespace cath { namespace rslv { struct html_hit; } }

namespace cath {
	namespace rslv {

		/// \brief Type alias for an optional crh_segment_spec
		using crh_segment_spec_opt = boost::optional<crh_segment_spec>;

		/// \brief Type alias for a vector of alnd_rgn
		using alnd_rgn_vec = std::vector<alnd_rgn>;

		/// \brief Type alias for an optional alnd_rgn_vec
		using alnd_rgn_vec_opt = boost::optional<alnd_rgn_vec>;

		/// \brief Type alias for a vector of hit_seg objects
		using hit_seg_vec = std::vector<hit_seg>;

		/// \brief Type alias for an optional hit_seg
		using hit_seg_opt = boost::optional<hit_seg>;

		/// \brief Type alias for the type to be used to index hits
		using hitidx_t = unsigned int;

		/// \brief Type alias for a vector of hitidx_t values
		using hitidx_vec = std::vector<hitidx_t>;

		/// \brief Type alias for a vector of calc_hit objects
		using calc_hit_vec = std::vector<calc_hit>;

		/// \brief Type alias for a vector of full_hit objects
		using full_hit_vec = std::vector<full_hit>;

		/// \brief Type alias for a vector of res_arrow objects
		using res_arrow_vec = std::vector<res_arrow>;

		/// \brief Type alias for a pair of res_arrows
		using res_arr_res_arr_pair = std::pair<res_arrow, res_arrow>;

		/// \brief Type alias for a pair of res_arrows
		using res_arr_res_arr_pair_vec = std::vector<res_arr_res_arr_pair>;

		/// \brief Type alias for an optional res_arrow
		using res_arrow_opt = boost::optional<res_arrow>;

		/// \brief Type alias for a pair of res_arrow_opt values
		using seg_boundary_pair     = std::pair<res_arrow_opt, res_arrow_opt>;

		/// \brief Type alias for a vector of seg_boundary_pair values
		using seg_boundary_pair_vec = std::vector<seg_boundary_pair>;

		/// \brief Type alias for an optional seg_boundary_pair_vec
		using seg_boundary_pair_vec_opt = boost::optional<seg_boundary_pair_vec>;

		/// \brief Type alias for the type to be used to index residues
		using residx_t = unsigned int;

		/// \brief Type alias for an optional residx_t
		using residx_opt = boost::optional<residx_t>;

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

		/// \brief Type alias for an optional resscr_t
		using resscr_opt = boost::optional<resscr_t>;

		/// \brief Type alias for an optional scored_arch_proxy object
		using scored_arch_proxy_opt = boost::optional<scored_arch_proxy>;

		/// \brief Type alias for a vector of scored_arch_proxy objects
		using scored_arch_proxy_vec = std::vector<scored_arch_proxy>;

		/// \brief The initial score before any hits have been added
		constexpr resscr_t INIT_SCORE = 0.0;

		/// \brief Type alias for a vector of html_hits
		using html_hit_vec = std::vector<html_hit>;

	} // namespace rslv
} // namespace cath

#endif
