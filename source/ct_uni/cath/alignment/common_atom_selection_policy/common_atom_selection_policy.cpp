/// \file
/// \brief The common_atom_selection_policy class definitions

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

#include "common_atom_selection_policy.hpp"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/lexical_cast.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/common_atom_selection_policy/common_atom_select_ca_policy.hpp"
#include "cath/alignment/common_atom_selection_policy/common_atom_select_cb_policy.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/coord_list.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::geom;
using namespace ::std;

using ::boost::assign::ptr_push_back;
using ::boost::ptr_vector;
using ::std::make_unique;

/// \brief NVI pass-through method to the concrete class's do_select_common_atoms() which defines how the policy extracts common atoms from a pair of residues
void common_atom_selection_policy::append_common_atoms_to_coord_lists(coord_list_coord_list_pair &prm_coord_lists, ///< The previous common coord_lists to which the new selections should be appended
                                                                      const residue              &prm_residue_a,   ///< The first  residue from which common atoms should be extracted
                                                                      const residue              &prm_residue_b    ///< The second residue from which common atoms should be extracted
                                                                      ) const {
	// Sanity check the inputs
	if ( prm_coord_lists.first.size() != prm_coord_lists.second.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The common coord_lists passed to append_common_atoms_to_coord_lists() must be of equal length"));
	}

	do_append_common_atoms_to_coord_lists( prm_coord_lists, prm_residue_a, prm_residue_b );

	// Sanity check the altered coord_lists
	if ( prm_coord_lists.first.size() != prm_coord_lists.second.size() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("common_atom_selection_policy has appended an unequal number of common coordinates to the two lists"));
	}
}

/// \brief NVI pass-through method to the concrete class's do_get_descriptive_name() which defines a string that descriptively names the policy
string common_atom_selection_policy::get_descriptive_name() const {
	return do_get_descriptive_name();
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<common_atom_selection_policy> common_atom_selection_policy::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief An NVI pass-through to the concrete class's do_less_than_with_same_dynamic_type(),
///        which defines the less-than operator when the argument's known to have the same dynamic type
bool common_atom_selection_policy::less_than_with_same_dynamic_type(const common_atom_selection_policy &prm_common_atom_selection_policy ///< TODOCUMENT
                                                     ) const {
	assert( typeid( *this ) == typeid( prm_common_atom_selection_policy ) );
	return do_less_than_with_same_dynamic_type( prm_common_atom_selection_policy );
}

/// \brief Convenience function for creating a new pair of coord_lists that the common_atom_selection_policy selects from the two residues
///
/// \relates common_atom_selection_policy
coord_list_coord_list_pair cath::align::select_common_atoms(const common_atom_selection_policy &prm_comm_atom_seln_pol, ///< The common_atom_selection_policy to make the selections
                                                            const residue                      &prm_residue_a,          ///< The first  residue from which common atoms should be extracted
                                                            const residue                      &prm_residue_b           ///< The second residue from which common atoms should be extracted
                                                            ) {
	// Create a new pair of coord_lists, perform the selections and then return the results
	pair<coord_list, coord_list> new_selections;
	prm_comm_atom_seln_pol.append_common_atoms_to_coord_lists( new_selections, prm_residue_a, prm_residue_b );
	return new_selections;
}

/// \brief Factory function that generates a list of all possible different common_atom_selection_policy objects
///
/// \relates common_atom_selection_policy
ptr_vector<common_atom_selection_policy> cath::align::get_all_common_atom_selection_policies() {
	ptr_vector<common_atom_selection_policy> all_policies;
	ptr_push_back<common_atom_select_ca_policy>( all_policies )();
	ptr_push_back<common_atom_select_cb_policy>( all_policies )();
	return all_policies;
}

/// \brief Factory function that creates the default common_atom_selection_policy (which selects alpha carbons)
///
/// \relates common_atom_selection_policy
unique_ptr<common_atom_selection_policy> cath::align::make_default_common_atom_selection_policy() {
	return { make_unique<common_atom_select_ca_policy>() };
}

/// \brief Return whether the specified common_atom_selection_policy is of the default type
///
/// ATM, this compares the get_descriptive_name, though it could equally use a < operator if one is added to the
/// common_atom_selection_policy hierarchy
///
/// \relates common_atom_selection_policy
bool cath::align::is_default_policy(const common_atom_selection_policy &prm_comm_atom_seln_pol ///< The prm_comm_atom_seln_pol to be tested for default-ness
                                    ) {
	return ( make_default_common_atom_selection_policy()->get_descriptive_name() == prm_comm_atom_seln_pol.get_descriptive_name() );
}
