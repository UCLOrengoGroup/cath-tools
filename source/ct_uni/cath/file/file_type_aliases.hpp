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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_FILE_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_UNI_FILE_FILE_TYPE_ALIASES_HPP

#include <boost/filesystem/path.hpp>

#include "cath/common/type_aliases.hpp"
#include "cath/file/pdb/pdb_atom_parse_status.hpp"

#include <utility>
#include <vector>

namespace cath { class amino_acid; }
namespace cath { class chain_label; }
namespace cath { class residue_id; }
namespace cath { namespace file { class backbone_complete_indices; } }
namespace cath { namespace file { class hmmer_scores_entry; } }
namespace cath { namespace file { class name_set; } }
namespace cath { namespace file { class name_set_list; } }
namespace cath { namespace file { class pdb; } }
namespace cath { namespace file { class pdb_atom; } }
namespace cath { namespace file { class pdb_list; } }
namespace cath { namespace file { class pdb_residue; } }
namespace cath { namespace file { class prc_scores_entry; } }
namespace cath { namespace file { class sec_file_record; } }
namespace cath { namespace file { class ssap_scores_entry; } }
namespace cath { namespace file { enum class data_file : unsigned int; } }
namespace cath { namespace file { enum class residue_makeup : bool; } }

namespace cath {
	namespace file {

		/// \brief Type alias for a vector of backbone_complete_indices values
		using backbone_complete_indices_vec = std::vector<backbone_complete_indices>;

		/// \brief Type alias for a vector of name_set
		using name_set_vec = std::vector<name_set>;

		/// \brief Type alias for a vector of residue_makeup values
		using residue_makeup_vec = std::vector<residue_makeup>;

		/// \brief TODOCUMENT
		using pdb_atom_vec = std::vector<pdb_atom>;

		/// \brief TODOCUMENT
		using pdb_residue_vec = std::vector<pdb_residue>;

		/// \brief Type alias for a pair of pdb_list and name_set_list
		using pdb_list_name_set_list_pair = std::pair<pdb_list, name_set_list>;

		/// \brief TODOCUMENT
		using pdb_vec = std::vector<pdb>;

		/// \brief Type alias for a pair of pdb and size_vec
		using pdb_size_vec_pair = std::pair<pdb, size_vec>;


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


		/// \brief Type alias for pair of residue_id and pdb_atom
		using res_id_pdb_atom_pair     = std::pair<residue_id, pdb_atom>;

		/// \brief Type alias for a vector of res_id_pdb_atom_pair values
		using res_id_pdb_atom_pair_vec = std::vector<res_id_pdb_atom_pair>;


		/// \brief Type alias for a reference_wrapper of a const ssap_scores_entry>
		using ssap_scores_entry_cref     = std::reference_wrapper<const ssap_scores_entry>;

		/// \brief Type alias for an optional reference_wrapper of a const ssap_scores_entry>
		using ssap_scores_entry_cref_opt = boost::optional<ssap_scores_entry_cref>;

		/// \brief Type alias for a pair of residue_id and pdb_atom
		using resid_atom_pair = std::pair<residue_id, pdb_atom>;

		/// \brief Type alias for a tuple of pdb_atom_parse_status, residuestring_name and amino_acid
		using status_string_aa_tuple = std::tuple<pdb_atom_parse_status, std::string, amino_acid>;


		/// \brief Type alias for a vector of sec_file_record objects
		using sec_file_record_vec = std::vector<sec_file_record>;

	} // namespace file
} // namespace cath

#endif
