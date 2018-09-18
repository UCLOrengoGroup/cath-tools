/// \file
/// \brief The tm_score class definitions

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

#include "tm_score.hpp"

#include <boost/logic/tribool.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/join.hpp>
#include <boost/range/numeric.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
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
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::score::detail;
using namespace cath::sup;
using namespace std;

using boost::inner_product;
using boost::numeric_cast;
using boost::tribool;

//BOOST_CLASS_EXPORT(tm_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> tm_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool tm_score::do_higher_is_better() const {
	return true;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value tm_score::do_calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                   const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                   const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                   ) const {
	// Extract the common coordinates to be chosen
	const pair<coord_list_vec, coord_list_vec> common_coords = the_coord_handler.get_common_coords_by_residue(
		prm_alignment,
		prm_protein_a,
		prm_protein_b
	);

	return max(
		score_for_target_length( common_coords, numeric_cast<score_value>( prm_protein_a.get_length() ) ),
		score_for_target_length( common_coords, numeric_cast<score_value>( prm_protein_b.get_length() ) )
	);
}

/// \brief Concrete implementation that describes what this score means
string tm_score::do_description() const {
	return "DESCRIPTION HERE PLEASE";
}

/// \brief TODOCUMENT
string tm_score::do_id_name() const {
	return "TM-score";
}

/// \brief TODOCUMENT
str_bool_pair_vec tm_score::do_short_name_suffixes() const {
	return the_coord_handler.short_name_suffixes();
}

/// \brief Concrete implementation providing long name
string tm_score::do_long_name() const {
	return "!!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!! !!LONG_NAME_HERE_PLEASE!!";
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string tm_score::do_reference() const {
	return { "Zhang Y, Skolnick J. Scoring function for automated assessment of protein structure template quality. Proteins. 2004 Dec 1;57(4):702-10. PubMed PMID: 15476259." };
}

/// \brief TODOCUMENT
bool tm_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                   ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief TODOCUMENT
score_value tm_score::score_for_target_length(const pair<coord_list_vec, coord_list_vec> &prm_common_coords_by_residue,  ///< TODOCUMENT
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

/// \brief Ctor for tm_score that allows the caller to specify the protein_only_length_getter, common_residue_selection_policy and common_atom_selection_policy
tm_score::tm_score(const common_residue_selection_policy &prm_comm_res_seln_pol, ///< The policy to use for selecting common residues
                   const common_atom_selection_policy    &prm_comm_atom_seln_pol ///< The policy to use for selecting common atoms
                   ) : the_coord_handler ( prm_comm_res_seln_pol, prm_comm_atom_seln_pol  ) {
}

/// \brief TODOCUMENT
const score_common_coord_handler & tm_score::get_score_common_coord_handler() const {
	return the_coord_handler;
}

/// \brief TODOCUMENT
///
/// \relates tm_score
bool cath::score::operator<(const tm_score &prm_tm_score_a, ///< TODOCUMENT
                            const tm_score &prm_tm_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( prm_tm_score_a, prm_tm_score_b );
	the_helper.register_comparison_field( &tm_score::get_score_common_coord_handler );
	return final_less_than_result( the_helper );
}
