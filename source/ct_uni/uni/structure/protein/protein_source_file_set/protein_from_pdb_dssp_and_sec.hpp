/// \file
/// \brief The protein_from_pdb_dssp_and_sec class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_FROM_PDB_DSSP_AND_SEC_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_PROTEIN_PROTEIN_SOURCE_FILE_SET_PROTEIN_FROM_PDB_DSSP_AND_SEC_HPP

#include "structure/protein/protein_source_file_set/protein_source_file_set.hpp"
#include "file/pdb/dssp_skip_policy.hpp"

namespace cath {

/// \brief Concrete protein_source_file_set for reading each protein from a PDB file, a DDSP file and a sec file
	class protein_from_pdb_dssp_and_sec final : public protein_source_file_set {
	private:
		/// \brief Whether to limit the protein to residues that are found in the DSSP file
		file::dssp_skip_policy the_dssp_skip_policy;

		std::unique_ptr<protein_source_file_set> do_clone() const final;

		file::data_file_vec do_get_file_set() const final;

		file::data_file do_get_primary_file() const final;

		protein_file_combn do_get_protein_file_combn() const final;

		bool do_makes_ssap_ready_protein() const final;

		protein do_read_files(const file::data_file_path_map &,
		                      const std::string &,
		                      std::ostream &) const final;

	public:
		explicit protein_from_pdb_dssp_and_sec(const file::dssp_skip_policy & = file::dssp_skip_policy::DONT_SKIP__DONT_BREAK_ANGLES);
	};

} // namespace cath

#endif
