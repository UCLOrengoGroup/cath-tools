/// \file
/// \brief The rmsd_score class definitions

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

#include "rmsd_score.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/serialization/export.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/less_than_helper.hpp"
#include "exception/out_of_range_exception.hpp"
#include "structure/geometry/coord_list.hpp"
#include "superposition/superposition.hpp"

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::score::detail;
using namespace cath::sup;
using namespace std;

BOOST_CLASS_EXPORT(rmsd_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> rmsd_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower RMSD score is better
tribool rmsd_score::do_higher_is_better() const {
	return false;
}

/// \brief Concrete implementation for calculating the RMSD of an alignment
score_value rmsd_score::do_calculate(const alignment &arg_alignment, ///< The pair alignment to be scored
                                     const protein   &arg_protein_a, ///< The protein associated with the first  half of the alignment
                                     const protein   &arg_protein_b  ///< The protein associated with the second half of the alignment
                                     ) const {
	// Extract the common coordinates to be chosen
	const pair<coord_list, coord_list> common_coords = the_coord_handler.get_common_coords(
		arg_alignment,
		arg_protein_a,
		arg_protein_b
	);

	// Check that there are some coords
	const size_t num_common_coords = common_coords.first.size();
	if ( num_common_coords == 0 ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"Cannot calculate "
			+ do_id_name()
			+ " for alignment from which no common coords are extracted"
		));
	}

	// Calculate and return the RMSD
	return calc_pairwise_superposition_rmsd(
		common_coords.first,
		common_coords.second
	);
}

/// \brief Concrete implementation that describes what this score means
string rmsd_score::do_description() const {
	return "The root of the mean squared deviation between equivalent atoms"
	       + the_coord_handler.description_brackets_string()
	       + ". This is measured in angstroms.";
}

/// \brief TODOCUMENT
string rmsd_score::do_id_name() const {
	return "RMSD";
}

/// \brief TODOCUMENT
str_bool_pair_vec rmsd_score::do_short_name_suffixes() const {
	return the_coord_handler.short_name_suffixes();
}

/// \brief Concrete implementation providing long name
string rmsd_score::do_long_name() const {
	return "Root Mean Squared Deviation" + the_coord_handler.long_suffix_string();
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> rmsd_score::do_build_from_short_name_spec(const string &arg_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                         ) const {
//	cerr << "Should build a rmsd_score from string \"" << arg_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool rmsd_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &arg_aligned_pair_score ///< TODOCUMENT
                                                     ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( arg_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Ctor for rmsd_score that allows the caller to specify the common_residue_selection_policy and common_atom_selection_policy
rmsd_score::rmsd_score(const common_residue_selection_policy &arg_comm_res_seln_pol, ///< The policy to use for selecting common residues
                       const common_atom_selection_policy    &arg_comm_atom_seln_pol   ///< The policy to use for selecting common atoms
                       ) : the_coord_handler(arg_comm_res_seln_pol, arg_comm_atom_seln_pol) {

}

/// \brief Pass-through method to provide public access to the score_the_coord_handler's short_suffix_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string rmsd_score::short_suffix_string() const {
	return the_coord_handler.short_suffix_string();
}

/// \brief Pass-through method to provide public access to the score_the_coord_handler's long_suffix_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string rmsd_score::long_suffix_string() const {
	return the_coord_handler.long_suffix_string();
}

/// \brief Pass-through method to provide public access to the score_the_coord_handler's description_brackets_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string rmsd_score::description_brackets_string() const {
	return the_coord_handler.description_brackets_string();
}

/// \brief TODOCUMENT
const score_common_coord_handler & rmsd_score::get_score_common_coord_handler() const {
	return the_coord_handler;
}

/// \brief TODOCUMENT
///
/// \relates rmsd_score
bool cath::score::operator<(const rmsd_score &arg_rmsd_score_a, ///< TODOCUMENT
                            const rmsd_score &arg_rmsd_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( arg_rmsd_score_a, arg_rmsd_score_b );
	the_helper.register_comparison_field( &rmsd_score::get_score_common_coord_handler );
	return final_less_than_result( the_helper );
}
