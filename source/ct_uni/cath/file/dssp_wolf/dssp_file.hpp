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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_DSSP_FILE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_DSSP_FILE_HPP

#include <optional>
#include <string>
#include <vector>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/name_set/name_set.hpp"
#include "cath/file/pdb/dssp_skip_policy.hpp"
#include "cath/structure/structure_type_aliases.hpp"

// clang-format off
namespace cath { class protein; }
namespace cath { class residue; }
namespace cath::file { class pdb; }
// clang-format on

namespace cath::file {

	/// \brief Represent the data parsed from a DSSP file
	class dssp_file final {
	private:
		residue_vec dssp_residues;

	public:
		using const_iterator = residue_vec_citr;

		explicit dssp_file(residue_vec);

		[[nodiscard]] size_t         get_num_residues() const;
		[[nodiscard]] const residue &get_residue_of_index( const size_t & ) const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	protein protein_from_dssp_and_pdb(const dssp_file &,
	                                  const pdb &,
	                                  const dssp_skip_policy & = dssp_skip_policy::DONT_SKIP__DONT_BREAK_ANGLES,
	                                  const file::name_set & = file::name_set{},
	                                  const ostream_ref_opt & = ::std::nullopt );

	residue_id_vec get_residue_ids(const dssp_file &,
	                               const bool &);

	size_t get_num_non_null_residues(const dssp_file &);

} // namespace cath::file

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_DSSP_WOLF_DSSP_FILE_HPP
