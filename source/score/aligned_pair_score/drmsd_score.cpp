/// \file
/// \brief The drmsd_score class definitions

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

#include "drmsd_score.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/serialization/export.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "exception/out_of_range_exception.hpp"
#include "structure/geometry/coord_list.hpp"
#include "superposition/superposition.hpp"

#include <cmath>

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::score::detail;
using namespace std;

using boost::numeric_cast;

BOOST_CLASS_EXPORT(drmsd_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> drmsd_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower dRMSD score is better
tribool drmsd_score::do_higher_is_better() const {
	return false;
}

/// \brief Concrete implementation for calculating the dRMSD of an alignment
score_value drmsd_score::do_calculate(const alignment &arg_alignment, ///< The pair alignment to be scored
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

	// Loop over the pairs
	double total_squared_deviation = 0.0;
	for (const size_t &coord_ctr_1 : indices( num_common_coords ) ) {
		for (size_t coord_ctr_2 = coord_ctr_1 + 1; coord_ctr_2 < num_common_coords; ++coord_ctr_2) {

			// Get the two equivalent distances between these pairs of atoms and add
			// the square of the difference to total_squared_deviation
			const double distance_a  = distance_between_points( common_coords.first[  coord_ctr_1 ], common_coords.first[  coord_ctr_2 ] );
			const double distance_b  = distance_between_points( common_coords.second[ coord_ctr_1 ], common_coords.second[ coord_ctr_2 ] );
			total_squared_deviation += (distance_a - distance_b) * (distance_a - distance_b);
		}
	}

	// Return the square root of the mean squared deviation
	const size_t num_values             = num_common_coords * (num_common_coords - 1) / 2;
	const double mean_squared_deviation = total_squared_deviation / numeric_cast<double>( num_values );
	return sqrt( mean_squared_deviation );
}

/// \brief Concrete implementation that describes what this score means
string drmsd_score::do_description() const {
	return "The root of the mean squared deviation between equivalent atoms"
	       + the_coord_handler.description_brackets_string()
	       + ". This is measured in angstroms.";
}

/// \brief TODOCUMENT
string drmsd_score::do_id_name() const {
	return "dRMSD";
}

/// \brief TODOCUMENT
str_bool_pair_vec drmsd_score::do_short_name_suffixes() const {
	return the_coord_handler.short_name_suffixes();
}

/// \brief Concrete implementation providing long name
string drmsd_score::do_long_name() const {
	return "Distance Root Mean Squared Deviation" + the_coord_handler.long_suffix_string();
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> drmsd_score::do_build_from_short_name_spec(const string &arg_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                          ) const {
//	cerr << "Should build a drmsd_score from string \"" << arg_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool drmsd_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &arg_aligned_pair_score ///< TODOCUMENT
                                                      ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( arg_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Ctor for drmsd_score that allows the caller to specify the common_residue_selection_policy and common_atom_selection_policy
drmsd_score::drmsd_score(const common_residue_selection_policy &arg_comm_res_seln_pol, ///< The policy to use for selecting common residues
                         const common_atom_selection_policy    &arg_comm_atom_seln_pol ///< The policy to use for selecting common atoms
                         ) : the_coord_handler( arg_comm_res_seln_pol, arg_comm_atom_seln_pol ) {
}

/// \brief Pass-through method to provide public access to the score_common_coord_handler's short_suffix_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string drmsd_score::short_suffix_string() const {
	return the_coord_handler.short_suffix_string();
}

/// \brief Pass-through method to provide public access to the score_common_coord_handler's long_suffix_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string drmsd_score::long_suffix_string() const {
	return the_coord_handler.long_suffix_string();
}

/// \brief Pass-through method to provide public access to the score_common_coord_handler's description_brackets_string() to help with
///        any other aligned_pair_scores that are implemented in terms of this class.
string drmsd_score::description_brackets_string() const {
	return the_coord_handler.description_brackets_string();
}

/// \brief TODOCUMENT
const score_common_coord_handler & drmsd_score::get_score_common_coord_handler() const {
	return the_coord_handler;
}

/// \brief TODOCUMENT
///
/// \relates drmsd_score
bool cath::score::operator<(const drmsd_score &arg_drmsd_score_a, ///< TODOCUMENT
                            const drmsd_score &arg_drmsd_score_b  ///< TODOCUMENT
                            ) {
	return arg_drmsd_score_a.get_score_common_coord_handler() < arg_drmsd_score_b.get_score_common_coord_handler();
}

