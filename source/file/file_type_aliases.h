/// \file
/// \brief The file type_aliases header

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

#ifndef FILE_TYPE_ALIASES_H_INCLUDED
#define FILE_TYPE_ALIASES_H_INCLUDED

#include "common/type_aliases.h"

#include <utility>
#include <vector>

namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class pdb_atom; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace file { class pdb_residue; } }

namespace cath {
	namespace file {
		
		/// \brief TODOCUMENT
		using pdb_atom_vec = std::vector<pdb_atom>;

		/// \brief TODOCUMENT
		using pdb_residue_vec = std::vector<pdb_residue>;

		/// \brief TODOCUMENT
		using pdb_list_str_vec_pair = std::pair<pdb_list, str_vec>;

		/// \brief TODOCUMENT
		using pdb_vec = std::vector<pdb>;
	}
}

#endif
