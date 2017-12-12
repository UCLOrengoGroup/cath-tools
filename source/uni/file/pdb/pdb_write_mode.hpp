/// \file
/// \brief The pdb_write_mode class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_PDB_PDB_WRITE_MODE_H
#define _CATH_TOOLS_SOURCE_UNI_FILE_PDB_PDB_WRITE_MODE_H

namespace cath {
	namespace file {

		/// \brief Whether a PDB to be written is the last/only part of the PDB or is one of several
		///        (eg where multiple PDBs are being faked as different chains from the same PDB)
		enum class pdb_write_mode : bool {
			ONLY_OR_LAST_PDB, ///< The last/only part of the PDB
			MORE_TO_FOLLOW    ///< To be followed by more parts in the PDB
		};

	} // namespace file
} // namespace cath

#endif
