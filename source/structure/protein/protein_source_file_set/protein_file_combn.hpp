/// \file
/// \brief The protein_file_combn header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_FILE_COMBN_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_FILE_COMBN_H

#include <iosfwd>
#include <memory>

namespace cath { class protein_source_file_set; }

namespace cath {
	class protein_source_file_set;

    /// \brief Represent the different sets of files from which proteins can be read
    ///
    /// This is a mirror of the protein_source_file_set class hierarchy that is useful when
    /// handling program_options
	enum class protein_file_combn : char {
		WOLF_SEC,     ///< Reading a protein from Wolf and sec files
		PDB,          ///< Reading a protein from PDB file
		PDB_DSSP_SEC, ///< Reading a protein from PDB, DSSP and sec files
		PDB_DSSP      ///< Reading a protein from PDB and DSSP (currently experimental)
	};

	std::unique_ptr<const protein_source_file_set> get_protein_source_file_set(const protein_file_combn &);

    std::istream & operator>>(std::istream &,
                              protein_file_combn &);

    std::ostream & operator<<(std::ostream &,
                              const protein_file_combn &);
} // namespace cath

#endif
