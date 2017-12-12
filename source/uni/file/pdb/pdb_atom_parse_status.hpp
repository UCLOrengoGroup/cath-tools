/// \file
/// \brief The pdb_atom_parse_status header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_PDB_PDB_ATOM_PARSE_STATUS_H
#define _CATH_TOOLS_SOURCE_UNI_FILE_PDB_PDB_ATOM_PARSE_STATUS_H

#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief The status arising from attempting to parse a pdb_atom from a string
		enum class pdb_atom_parse_status : char {
			OK,   ///< The record parses correctly
			SKIP, ///< The record fails with a minor error; skip this record but process further records.
			ABORT ///< The record fails with a major error; abort parse attempt.
		};

//		pdb_atom_parse_status str_to_pdb_rec(const std::string &);
//
//		std::istream & operator>>(std::istream &,
//		                          pdb_atom_parse_status &);
//
//		std::ostream & operator<<(std::ostream &,
//		                          const pdb_atom_parse_status &);

	} // namespace file
} // namespace cath

#endif
