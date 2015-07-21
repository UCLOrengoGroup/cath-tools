/// \file
/// \brief The common_atom_select_cb_policy class definitions

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

#include "common_atom_select_cb_policy.h"

#include "common/clone/make_uptr_clone.h"
#include "structure/geometry/coord_list.h"
#include "structure/protein/residue.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

/// \brief Concrete implementation for selecting common atoms by selecting the two residues' carbon beta atoms
void common_atom_select_cb_policy::do_append_common_atoms_to_coord_lists(coord_list_coord_list_pair &arg_coord_lists, ///< The previous common coord_lists to which the new selections should be appended
                                                                         const residue              &arg_residue_a,   ///< The first  residue from which common atoms should be extracted
                                                                         const residue              &arg_residue_b    ///< The second residue from which common atoms should be extracted
                                                                         ) const {
	arg_coord_lists.first.push_back ( arg_residue_a.get_carbon_beta_coord() );
	arg_coord_lists.second.push_back( arg_residue_b.get_carbon_beta_coord() );
}

/// Concrete implementation for providing a descriptive name of this policy
string common_atom_select_cb_policy::do_get_descriptive_name() const {
	return "cb_atoms";
}

/// \brief A standard do_clone method.
unique_ptr<common_atom_selection_policy> common_atom_select_cb_policy::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
bool common_atom_select_cb_policy::do_less_than_with_same_dynamic_type(const common_atom_selection_policy &/*arg_common_atom_selection_policy*/ ///< TODOCUMENT
                                                                       ) const {
	return false;
}
