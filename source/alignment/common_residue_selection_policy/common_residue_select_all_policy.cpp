/// \file
/// \brief The common_residue_select_all_policy class definitions

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

#include "common_residue_select_all_policy.h"

#include "common/clone/make_uptr_clone.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
size_vec common_residue_select_all_policy::do_select_common_residues(const alignment                    &/*arg_alignment*/, ///< TODOCUMENT
                                                                     const vector<alignment::size_type> &arg_indices,       ///< TODOCUMENT
                                                                     const alignment::size_type         &/*arg_entry_a*/,   ///< TODOCUMENT
                                                                     const alignment::size_type         &/*arg_entry_b*/    ///< TODOCUMENT
                                                                     ) const {
	size_vec index_indices;
	index_indices.reserve(arg_indices.size());
	for (size_t index_ctr = 0; index_ctr < arg_indices.size(); ++index_ctr) {
		index_indices.push_back(index_ctr);
	}
	return index_indices;
}

string common_residue_select_all_policy::do_get_descriptive_name() const {
	return "select_all";
}

/// \brief TODOCUMENT
unique_ptr<common_residue_selection_policy> common_residue_select_all_policy::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
bool common_residue_select_all_policy::do_less_than_with_same_dynamic_type(const common_residue_selection_policy &/*arg_common_residue_selection_policy*/ ///< TODOCUMENT
                                                                           ) const {
	return false;
}
