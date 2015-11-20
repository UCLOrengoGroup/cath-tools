/// \file
/// \brief The file type_aliases header

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

#ifndef FILE_TYPE_ALIASES_H_INCLUDED
#define FILE_TYPE_ALIASES_H_INCLUDED

#include "common/type_aliases.h"

#include <utility>
#include <vector>

namespace cath { namespace file { enum class data_file : unsigned int; } }
namespace cath { namespace file { class hmmer_scores_entry; } }
namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class pdb_atom; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace file { class pdb_residue; } }
namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { class ssap_scores_entry; } }

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


		/// \brief Type alias for vector of data_files
		using data_file_vec = std::vector<data_file>;

		/// \brief Type alias for set of data_files
		using data_file_set = std::set<data_file>;

		/// \brief Type alias for pair of data_file and path
		using data_file_path_pair = std::pair<data_file, boost::filesystem::path>;

		/// \brief Type alias for map from data_file to path
		using data_file_path_map = std::map<data_file, boost::filesystem::path>;


		/// \brief Type alias for a vector of ssap_scores_entry objects
		using ssap_scores_entry_vec = std::vector<ssap_scores_entry>;

		/// \brief Type alias for a vector of prc_scores_entry objects
		using prc_scores_entry_vec = std::vector<prc_scores_entry>;

		/// \brief Type alias for a vector of hmmer_scores_entry objects
		using hmmer_scores_entry_vec = std::vector<hmmer_scores_entry>;
	}
}

#endif
