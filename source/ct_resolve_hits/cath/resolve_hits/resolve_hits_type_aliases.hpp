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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE_HITS_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE_HITS_TYPE_ALIASES_HPP

#include <iosfwd>
#include <optional>
#include <vector>

#include "cath/common/type_aliases.hpp"
#include "cath/seq/seq_type_aliases.hpp"

// clang-format off
namespace cath::rslv { class calc_hit; }
namespace cath::rslv { class calc_hit_list; }
namespace cath::rslv { class crh_segment_spec; }
namespace cath::rslv { class full_hit; }
namespace cath::rslv { class full_hit_list; }
namespace cath::rslv { class scored_arch_proxy; }
namespace cath::rslv { class trim_spec; }
namespace cath::rslv { struct alnd_rgn; }
namespace cath::rslv { struct html_hit; }
namespace cath::rslv::detail { class full_hit_prune_builder; }
namespace cath::rslv::detail { class hits_processor; }
// clang-format on

namespace cath::rslv {

	/// \brief Type alias for a vector of alnd_rgn
	using alnd_rgn_vec                  = std::vector<alnd_rgn>;

	/// \brief Type alias for an optional alnd_rgn_vec
	using alnd_rgn_vec_opt              = ::std::optional<alnd_rgn_vec>;

	/// \brief Type alias for a vector of calc_hit objects
	using calc_hit_vec                  = std::vector<calc_hit>;

	/// \brief Type alias for calc_hit_vec's const_iterator type
	using calc_hit_vec_citr             = calc_hit_vec::const_iterator;

	/// \brief Type alias for a vector of calc_hit_vec_citr values
	using calc_hit_vec_citr_vec         = std::vector<calc_hit_vec_citr>;

	/// \brief Type alias for calc_hit_vec's iterator type
	using calc_hit_vec_itr              = calc_hit_vec::iterator;

	/// \brief Type alias for an optional crh_segment_spec
	using crh_segment_spec_opt          = ::std::optional<crh_segment_spec>;

	/// \brief Type alias for a vector of full_hit objects
	using full_hit_vec                  = std::vector<full_hit>;

	/// \brief Type alias for an optional full_hit_list
	using full_hit_list_opt             = ::std::optional<full_hit_list>;

	/// \brief Type alias for the type to be used to index hits
	using hitidx_t                      = unsigned int;

	/// \brief Type alias for a vector of hitidx_t values
	using hitidx_vec                    = std::vector<hitidx_t>;

	/// \brief Type alias for a vector of html_hits
	using html_hit_vec                  = std::vector<html_hit>;

	/// \brief Type alias for a pair of seq_arrow and hit index
	using res_arr_idx_pair              = std::pair<seq::seq_arrow, hitidx_t>;

	/// \brief Type alias for a vector of pairs of seq_arrow and hitidx_t
	using res_arr_idx_pair_vec          = std::vector<res_arr_idx_pair>;

	/// \brief Type alias for a pair of res_arrows
	using res_arr_res_arr_pair          = std::pair<seq::seq_arrow, seq::seq_arrow>;

	/// \brief Type alias for a pair of res_arrows
	using res_arr_res_arr_pair_vec      = std::vector<res_arr_res_arr_pair>;

	/// \brief Type alias for the type to be used for hits' scores
	using resscr_t                      = float;

	/// \brief Type alias for an optional resscr_t
	using resscr_opt                    = ::std::optional<resscr_t>;

	/// \brief Type alias for an optional scored_arch_proxy object
	using scored_arch_proxy_opt         = ::std::optional<scored_arch_proxy>;

	/// \brief Type alias for a vector of scored_arch_proxy objects
	using scored_arch_proxy_vec         = std::vector<scored_arch_proxy>;

	/// \brief Type alias for a pair of res_arrow_opt values
	using seg_boundary_pair             = std::pair<seq::res_arrow_opt, seq::res_arrow_opt>;

	/// \brief Type alias for a vector of seg_boundary_pair values
	using seg_boundary_pair_vec         = std::vector<seg_boundary_pair>;

	/// \brief Type alias for an optional seg_boundary_pair_vec
	using seg_boundary_pair_vec_opt     = ::std::optional<seg_boundary_pair_vec>;

	/// \brief Type alias for a pair of string and calc_hit_list
	using str_calc_hit_list_pair        = std::pair<std::string, calc_hit_list>;

	/// \brief Type alias for a vector of str_calc_hit_list_pair values
	using str_calc_hit_list_pair_vec    = std::vector<str_calc_hit_list_pair>;

	/// \brief Type alias for a pair of string and full_hit
	using str_full_hit_pair             = std::pair<std::string, full_hit>;

	/// \brief Type alias for an optional str_full_hit_pair
	using str_full_hit_pair_opt         = ::std::optional<str_full_hit_pair>;

	/// \brief Type alias for a pair of string and hit data (non-const) reference
	using str_hits_builder_ref_pair     = std::pair<std::string, detail::full_hit_prune_builder &>;

	/// \brief Type alias for an optional str_hits_ref_pair
	using str_hits_builder_ref_pair_opt = ::std::optional<str_hits_builder_ref_pair>;

	/// \brief Type alias for an optional trim_spec
	using trim_spec_opt                 = ::std::optional<trim_spec>;


	/// \brief The initial score before any hits have been added
	inline constexpr resscr_t INIT_SCORE = 0.0;

	namespace detail {

		/// \brief Type alias for a clone_ptr to a hits_processor
		using hits_processor_clptr = common::clone_ptr<hits_processor>;

		/// \brief Type alias for a vector of clone_ptrs to hits_processors
		using hits_processor_clptr_vec = common::clptr_vec<hits_processor>;

		/// \brief Type alias for a unique_ptr to a hits_processor
		using hits_processor_uptr = std::unique_ptr<hits_processor>;
	} // namespace detail

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_RESOLVE_HITS_TYPE_ALIASES_HPP
