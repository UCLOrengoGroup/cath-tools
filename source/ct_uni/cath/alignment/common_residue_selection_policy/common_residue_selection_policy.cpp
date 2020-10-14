/// \file
/// \brief The common_residue_selection_policy class definitions

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

#include "common_residue_selection_policy.hpp"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_all_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.hpp"
#include "cath/alignment/common_residue_selection_policy/common_residue_select_min_score_policy.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/common/cpp14/make_unique.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/size_t_literal.hpp"

using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::std;

using ::boost::adaptors::filtered;
using ::boost::assign::ptr_push_back;
using ::boost::irange;
using ::boost::lexical_cast;
using ::boost::ptr_vector;

/// \brief TODOCUMENT
vector<alignment::size_type> common_residue_selection_policy::select_common_residues(const alignment            &prm_alignment, ///< TODOCUMENT
                                                                                     const alignment::size_type &prm_entry_a,   ///< TODOCUMENT
                                                                                     const alignment::size_type &prm_entry_b    ///< TODOCUMENT
                                                                                     ) const {
	using aln_size_type = alignment::size_type;
	using aln_size_vec  = vector<aln_size_type>;

	// Sanity check both entries are valid
	const alignment::size_type num_entries = prm_alignment.num_entries();
	if (prm_entry_a >= num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Argument entry_a of " + lexical_cast<string>(prm_entry_a) +
			" is invalid for an argument of " + lexical_cast<string>(num_entries) + " entries"
		));
	}
	if (prm_entry_b >= num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Argument entry_b of " + lexical_cast<string>(prm_entry_b) +
			" is invalid for an argument of " + lexical_cast<string>(num_entries) + " entries"
		));
	}

	// Grab the indices of the positions that are in common between the entries prm_entry_a and prm_entry_b
	const auto original_indices = copy_build<aln_size_vec>(
		indices( prm_alignment.length() )
			| filtered(
				[&] (const size_t &x) {
					return (
						has_position_of_entry_of_index( prm_alignment, prm_entry_a, x )
						&&
						has_position_of_entry_of_index( prm_alignment, prm_entry_b, x )
					);
				}
			)
	);

	// Grab the results from the concrete class's implementation of the do_select_common_residues() method
	const size_vec common_coords_raw = do_select_common_residues(
		prm_alignment,
		original_indices,
		prm_entry_a,
		prm_entry_b
	);

	// Convert the indices back to indices that correspond to the original alignment
	aln_size_vec common_coords;
	common_coords.reserve(common_coords_raw.size());
	for (const size_t &prm_common_coord_raw : common_coords_raw) {
		if (prm_common_coord_raw >= original_indices.size()) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("do_select_common_residues produced raw common coord index out of range"));
		}
		common_coords.push_back(original_indices[prm_common_coord_raw]);
	}

	// Sanity check that the common_coords values are strictly increasing
	if ( ! common_coords.empty() ) {
		for (const size_t &common_residue_ctr : irange( 1_z, common_coords.size() ) ) {
			if (common_coords[common_residue_ctr] <= common_coords[common_residue_ctr-1]) {
				BOOST_THROW_EXCEPTION(runtime_error_exception("common_residue_selection_policy generated common coords that are not strictly increasing"));
			}
		}
	}

	return common_coords;
}

/// \brief TODOCUMENT
string common_residue_selection_policy::get_descriptive_name() const {
	return do_get_descriptive_name();
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<common_residue_selection_policy> common_residue_selection_policy::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief An NVI pass-through to the concrete class's do_less_than_with_same_dynamic_type(),
///        which defines the less-than operator when the argument's known to have the same dynamic type
bool common_residue_selection_policy::less_than_with_same_dynamic_type(const common_residue_selection_policy &prm_common_residue_selection_policy ///< TODOCUMENT
                                                                       ) const {
	assert( typeid( *this ) == typeid( prm_common_residue_selection_policy ) );
	return do_less_than_with_same_dynamic_type( prm_common_residue_selection_policy );
}

/// \brief Factory function that generates a list of all possible different common_residue_selection_policy objects
///
/// \relates common_residue_selection_policy
ptr_vector<common_residue_selection_policy> cath::align::get_all_common_residue_selection_policies() {
	ptr_vector<common_residue_selection_policy> all_policies;
	ptr_push_back< common_residue_select_all_policy                >( all_policies )(      );
	ptr_push_back< common_residue_select_best_score_percent_policy >( all_policies )(      );
	ptr_push_back< common_residue_select_min_score_policy          >( all_policies )( 0.01 );
	return all_policies;
}

/// \brief TODOCUMENT
///
/// \relates common_residue_selection_policy
vector<alignment::size_type> cath::align::select_common_residues_of_pair_alignment(const common_residue_selection_policy &prm_policy,   ///< TODOCUMENT
                                                                                   const alignment                       &prm_alignment ///< TODOCUMENT
                                                                                   ) {
	return prm_policy.select_common_residues(
		prm_alignment,
		alignment::PAIR_A_IDX,
		alignment::PAIR_B_IDX
	);
}


/// \brief Factory function that creates the default common_residue_selection_policy (which selects all aligned residue pairs)
///
/// \relates common_residue_selection_policy
unique_ptr<common_residue_selection_policy> cath::align::make_default_common_residue_selection_policy() {
	return { common::make_unique<common_residue_select_all_policy>() };
}

/// \brief Return whether the specified common_residue_selection_policy is of the default type
///
/// ATM, this compares the get_descriptive_name, though it could equally use a < operator if one is added to the
/// common_residue_selection_policy hierarchy
///
/// \relates common_residue_selection_policy
bool cath::align::is_default_policy(const common_residue_selection_policy &prm_comm_res_seln_pol ///< TODOCUMENT
                                    ) {
	return ( make_default_common_residue_selection_policy()->get_descriptive_name() == prm_comm_res_seln_pol.get_descriptive_name() );
}

/// \brief Simple insertion operator for common_residue_selection_policy
///
//// \relates common_residue_selection_policy
ostream & cath::align::operator<<(ostream                               &prm_os,                             ///< The ostream to which the common_residue_selection_policy should be output
                                  const common_residue_selection_policy &prm_common_residue_selection_policy ///< The protein_file_combn to output
                                  ) {
	prm_os << ( "common_residue_selection_policy[" + prm_common_residue_selection_policy.get_descriptive_name() + "]" );
	return prm_os;
}
