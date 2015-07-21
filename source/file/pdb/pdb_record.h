/// \file
/// \brief The pdb_record header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef PDB_RECORD_H_INCLUDED
#define PDB_RECORD_H_INCLUDED

#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		enum class pdb_record {
			ATOM,  ///< TODOCUMENT
			HETATM ///< TODOCUMENT
		};

		pdb_record str_to_pdb_rec(const std::string &);

		std::istream & operator>>(std::istream &,
		                          pdb_record &);

		std::ostream & operator<<(std::ostream &,
		                          const pdb_record &);

	}
}

#endif
