/// \file
/// \brief The dssp_file class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_DSSP_WOLF_DSSP_FILE_H
#define _CATH_TOOLS_SOURCE_FILE_DSSP_WOLF_DSSP_FILE_H

#include "common/type_aliases.h"
#include "structure/structure_type_aliases.h"

#include <string>
#include <vector>

namespace cath { namespace file { class pdb; } }
namespace cath { class protein; }
namespace cath { class residue; }

namespace cath {
	namespace file {

		/// \brief Represent the data parsed from a DSSP file
		class dssp_file final {
		private:
			residue_vec dssp_residues;

		public:
			using const_iterator = residue_vec_citr;

			dssp_file(const residue_vec &);

			size_t get_num_residues() const;
			const residue & get_residue_of_index(const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		protein protein_from_dssp_and_pdb(const dssp_file &,
		                                  const pdb &,
		                                  const bool &,
		                                  const std::string & = "");

		residue_name_vec get_residue_names(const dssp_file &,
		                                   const bool &);

		size_t get_num_non_null_residues(const dssp_file &);

	} // namespace file

} // namespace cath

#endif
