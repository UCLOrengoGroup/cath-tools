/// \file
/// \brief The common_residue_select_all_policy class definitions

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

#include "common_residue_select_all_policy.hpp"

#include "common/boost_addenda/range/indices.hpp"
#include "common/clone/make_uptr_clone.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
size_vec common_residue_select_all_policy::do_select_common_residues(const alignment                    &/*prm_alignment*/, ///< TODOCUMENT
                                                                     const vector<alignment::size_type> &prm_indices,       ///< TODOCUMENT
                                                                     const alignment::size_type         &/*prm_entry_a*/,   ///< TODOCUMENT
                                                                     const alignment::size_type         &/*prm_entry_b*/    ///< TODOCUMENT
                                                                     ) const {
	size_vec index_indices;
	index_indices.reserve(prm_indices.size());
	for (const size_t &index_ctr : indices( prm_indices.size() ) ) {
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
bool common_residue_select_all_policy::do_less_than_with_same_dynamic_type(const common_residue_selection_policy &/*prm_common_residue_selection_policy*/ ///< TODOCUMENT
                                                                           ) const {
	return false;
}
