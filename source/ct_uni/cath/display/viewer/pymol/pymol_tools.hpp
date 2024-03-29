/// \file
/// \brief The pymol_tools class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL_PYMOL_TOOLS_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL_PYMOL_TOOLS_HPP

#include <cstddef>
#include <optional>
#include <string>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath { class residue_id; }
// clang-format on

namespace cath {

	/// \brief TODOCUMENT
	struct pymol_tools final {
		pymol_tools() = delete;
		~pymol_tools() = delete;

		static double pymol_size(const size_t &,
		                         const double &,
		                         const size_t &,
		                         const double &,
		                         const size_t &);

		static std::string parse_residue_name_for_pymol(const residue_name &);

		static str_vec parse_residue_names_for_pymol(const residue_name_vec &);

		static std::string pymol_res_seln_str(const std::string &,
		                                      const residue_id_vec &,
		                                      const str_opt & = ::std::nullopt);

	};

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_DISPLAY_VIEWER_PYMOL_PYMOL_TOOLS_HPP
