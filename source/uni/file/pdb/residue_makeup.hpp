/// \file
/// \brief The residue_makeup header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_PDB_RESIDUE_MAKEUP_H
#define _CATH_TOOLS_SOURCE_UNI_FILE_PDB_RESIDUE_MAKEUP_H

namespace cath {
	namespace file {

		/// \brief Represent whether the records for a residue included non-proper amino-acids
		///        (say in HETATM records in a PDB file)
		enum class residue_makeup : bool {
			ALL_PROPER_AMINO_ACIDS,     ///< All records for this residue had proper amino-acids
			SOME_NON_PROPER_AMINO_ACIDS ///< Some (or all) records for this residue had non-proper amino-acids
		};

	} // namespace file
} // namespace cath

#endif
