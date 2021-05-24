/// \file
/// \brief The wolf_file class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_WOLF_FILE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_WOLF_FILE_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/structure/structure_type_aliases.hpp"

#include <vector>

// clang-format off
namespace cath::file { class pdb; }
namespace cath { class protein; }
namespace cath { class residue; }
// clang-format on

namespace cath::file {

	/// \brief Represent the data parsed from a WOLF file
	class wolf_file final {
		residue_vec wolf_residues;

	public:
		using size_type = residue_vec::size_type;

		explicit wolf_file(residue_vec);

		[[nodiscard]] size_type      get_num_residues() const;
		[[nodiscard]] const residue &get_residue_of_index( const size_type & ) const;
	};

	protein protein_from_wolf(const wolf_file &);

} // namespace cath::file

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_WOLF_FILE_HPP
