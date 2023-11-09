/// \file
/// \brief The seq type_aliases header

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

#ifndef CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_TYPE_ALIASES_HPP
#define CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_TYPE_ALIASES_HPP

#include <optional>
#include <vector>

// clang-format off
namespace cath::seq { class seq_arrow; }
namespace cath::seq { class seq_seg; }
namespace cath::seq { class seq_seg_run; }
// clang-format on

namespace cath::seq {

	/// \brief Type alias for an optional seq_arrow
	using res_arrow_opt             = ::std::optional<seq_arrow>;

	/// \brief Type alias for a vector of seq_arrow objects
	using res_arrow_vec             = std::vector<seq_arrow>;

	/// \brief Type alias for the type to be used for indexing residue boundaries
	///        (within a seq_arrow object)
	using resarw_t                  = unsigned int;

	/// \brief Type alias for the type to be used to index residues
	using residx_t                  = resarw_t;

	/// \brief Type alias for an optional residx_t
	using residx_opt                = ::std::optional<residx_t>;

	/// \brief Type alias for a pair of residue indices
	using residx_residx_pair        = std::pair<residx_t, residx_t>;

	/// \brief Type alias for a vector of pairs of residue indices
	using residx_residx_pair_vec    = std::vector<residx_residx_pair>;

	/// \brief Type alias for a vector residx_t
	using residx_vec                = std::vector<residx_t>;

	/// \brief Type alias for a pair of seq_arrows
	using seq_arr_seq_arr_pair      = std::pair<seq_arrow, seq_arrow>;

	/// \brief Type alias for an optional seq_seg
	using seq_seg_opt               = ::std::optional<seq_seg>;

	/// \brief Type alias for a vector of seq_seg_opt
	using seq_seg_opt_vec           = std::vector<seq_seg_opt>;

	/// \brief Type alias for an optional seq_seg_run
	using seq_seg_run_opt           = ::std::optional<seq_seg_run>;

	/// \brief Type alias for a vector of seq_seg_run
	using seq_seg_run_vec           = std::vector<seq_seg_run>;

	/// \brief Type alias for a vector of seq_seg objects
	using seq_seg_vec               = std::vector<seq_seg>;

} // namespace cath::seq

#endif // CATH_TOOLS_SOURCE_CT_SEQ_CATH_SEQ_SEQ_TYPE_ALIASES_HPP
