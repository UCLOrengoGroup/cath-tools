/// \file
/// \brief The biocore type_aliases header

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

#ifndef _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_BIOCORE_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_BIOCORE_TYPE_ALIASES_HPP

#include <map>
#include <optional>
#include <set>
#include <vector>

// clang-format off
namespace cath { class chain_label; }
namespace cath { class residue_id; }
namespace cath { class residue_name; }
// clang-format on

namespace cath {

	/// \brief Type alias for a vector of chain_label objects
	using chain_label_vec                 = ::std::vector<chain_label>;

	/// \brief Type alias for a set of chain_label objects
	using chain_label_set                 = ::std::set<chain_label>;

	/// \brief TODOCUMENT
	using chain_label_opt                 = ::std::optional<chain_label>;

	/// \brief TODOCUMENT
	using residue_name_opt                = ::std::optional<residue_name>;

	/// \brief TODOCUMENT
	using residue_name_set                = ::std::set<residue_name>;

	/// \brief TODOCUMENT
	using residue_name_vec                = ::std::vector<residue_name>;

	/// \brief TODOCUMENT
	using residue_name_vec_itr            = residue_name_vec::iterator;

	/// \brief TODOCUMENT
	using residue_name_vec_citr           = residue_name_vec::const_iterator;

	/// \brief TODOCUMENT
	using residue_name_vec_vec            = ::std::vector<residue_name_vec>;

	/// \brief Type alias for a vector of residue_ids
	using residue_id_vec                  = ::std::vector<residue_id>;

	/// \brief Type alias for a vector of residue_id_vec
	using residue_id_vec_vec              = ::std::vector<residue_id_vec>;

	/// \brief Type alias for a map from chain_label keys to residue_id_vec values
	using chain_label_residue_id_vec_map  = ::std::map<chain_label, residue_id_vec>;

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_BIOCORE_TYPE_ALIASES_HPP
