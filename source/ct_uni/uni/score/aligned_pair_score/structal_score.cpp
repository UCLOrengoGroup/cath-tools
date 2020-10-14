/// \file
/// \brief The structal_score class definitions

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

#include "structal_score.hpp"

#include <boost/logic/tribool.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/join.hpp>
#include <boost/range/numeric.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "alignment/gap/alignment_gap.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/less_than_helper.hpp"
#include "score/length_getter/length_of_shorter_getter.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/structure_type_aliases.hpp"
#include "superposition/superposition.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::align::gap;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::score::detail;
using namespace cath::sup;
using namespace std;

using boost::inner_product;
using boost::numeric_cast;
using boost::range::join;
using boost::tribool;

//BOOST_CLASS_EXPORT(structal_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> structal_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool structal_score::do_higher_is_better() const {
	return true;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value structal_score::do_calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                         const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                         const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                         ) const {
	constexpr double numerator        = 20.00;
	constexpr double dest_denominator =  2.24;
	constexpr double gap_weight       = 10.00;

	// Extract the common coordinates to be chosen
	const pair<coord_list_vec, coord_list_vec> common_coords = the_coord_handler.get_common_coords_by_residue(
		prm_alignment,
		prm_protein_a,
		prm_protein_b
	);

	const score_value num_gaps = numeric_cast<score_value>( get_naive_num_gaps( prm_alignment ) );

	const auto superposed_second_coords = superpose_copy_second_coords_to_first(
		common_coords.first,
		common_coords.second
	);

	return inner_product(
		common_coords.first,
		superposed_second_coords,
		numeric_cast<score_value>( 0.0 ),
		plus<score_value>(),
		[&] (const coord_list &x, const coord_list &y) {
			const score_value mean_distance = calc_mean_deviation( x, y );
			const score_value fraction      = mean_distance / dest_denominator;
			const score_value final         = numerator / ( 1.0 + fraction * fraction );
			return final;
		}
	) - ( gap_weight * num_gaps );
}

/// \brief Concrete implementation that describes what this score means
string structal_score::do_description() const {
	return "DESCRIPTION HERE PLEASE";
}

/// \brief TODOCUMENT
string structal_score::do_id_name() const {
	return "structal";
}

/// \brief TODOCUMENT
str_bool_pair_vec structal_score::do_short_name_suffixes() const {
	return the_coord_handler.short_name_suffixes();
}

/// \brief Concrete implementation providing long name
string structal_score::do_long_name() const {
	return "!!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!!";
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string structal_score::do_reference() const {
	return { "!!STRUCTAL REFERENCE HERE PLEASE !! !!STRUCTAL REFERENCE HERE PLEASE !! !!STRUCTAL REFERENCE HERE PLEASE !! !!STRUCTAL REFERENCE HERE PLEASE !!" };
}

/// \brief TODOCUMENT
bool structal_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                         ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief TODOCUMENT
score_value structal_score::score_for_target_length(const pair<coord_list_vec, coord_list_vec> &prm_common_coords_by_residue,  ///< TODOCUMENT
                                                    const score_value                          &prm_target_length              ///< TODOCUMENT
                                                    ) {
	const score_value denominator = ( 1.24 * cbrt( prm_target_length - 15 ) ) - 1.8;
//	cerr << "For length " << prm_target_length << ", denominator is : " << denominator << endl;

	const auto superposed_second_coords = superpose_copy_second_coords_to_first(
		prm_common_coords_by_residue.first,
		prm_common_coords_by_residue.second
	);

	return inner_product(
		prm_common_coords_by_residue.first,
		superposed_second_coords,
		numeric_cast<score_value>( 0.0 ),
		plus<score_value>(),
		[&] (const coord_list &x, const coord_list &y) {

			const score_value mean_distance = calc_mean_deviation( x, y );
			const score_value fraction      = mean_distance / denominator;
			const score_value final         = 1.0 / ( 1.0 + fraction * fraction );
//			cerr << "mean_distance is : " << mean_distance << endl;
//			cerr << "fraction      is : " << fraction      << endl;
//			cerr << "final         is : " << final         << endl;
			return final;
		}
	) / prm_target_length;
}

/// \brief Ctor for structal_score that allows the caller to specify the protein_only_length_getter, common_residue_selection_policy and common_atom_selection_policy
structal_score::structal_score(const common_residue_selection_policy &prm_comm_res_seln_pol, ///< The policy to use for selecting common residues
                               const common_atom_selection_policy    &prm_comm_atom_seln_pol ///< The policy to use for selecting common atoms
                               ) : the_coord_handler ( prm_comm_res_seln_pol, prm_comm_atom_seln_pol  ) {
}

/// \brief TODOCUMENT
const score_common_coord_handler & structal_score::get_score_common_coord_handler() const {
	return the_coord_handler;
}

/// \brief TODOCUMENT
///
/// \relates structal_score
bool cath::score::operator<(const structal_score &prm_structal_score_a, ///< TODOCUMENT
                            const structal_score &prm_structal_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( prm_structal_score_a, prm_structal_score_b );
	the_helper.register_comparison_field( &structal_score::get_score_common_coord_handler );
	return final_less_than_result( the_helper );
}
