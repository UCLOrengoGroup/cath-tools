/// \file
/// \brief The protein_info class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PROTEIN_INFO_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PROTEIN_INFO_HPP

#include "cath/file/pdb/residue_makeup.hpp"

#include <vector>

namespace cath::file {

	/// \brief Represent extra information about a protein from the process of building it
	struct protein_info final {

		/// \brief The residue_makeup of each residue in the protein
		std::vector<residue_makeup> residue_makeups;
	};

} // namespace cath::file

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_PDB_PROTEIN_INFO_HPP
