/// \file
/// \brief The mi_score class definitions

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

#include "mi_score.hpp"

#include <boost/logic/tribool.hpp>
#include <boost/range/join.hpp>

#include "alignment/alignment.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/less_than_helper.hpp"
#include "score/length_getter/length_of_shorter_getter.hpp"

// #include <iostream> // ***** TEMPORARY *****

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::range::join;
using boost::tribool;

//BOOST_CLASS_EXPORT(mi_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> mi_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool mi_score::do_higher_is_better() const {
	return false;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value mi_score::do_calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                   const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                   const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                   ) const {
	const score_value w0                    = 1.5;
	const score_value one_plus_num_aligned  = 1 + get_length_score( num_aligned_getter,      prm_alignment, prm_protein_a, prm_protein_b );
	const score_value one_plus_length       = 1 + get_length_score( *full_length_getter_ptr,                prm_protein_a, prm_protein_b );
	const score_value one_plus_rmsd_over_w0 = 1 + rmsd.calculate( prm_alignment, prm_protein_a, prm_protein_b ) / w0;
	return 1 - (one_plus_num_aligned / ( one_plus_rmsd_over_w0 * one_plus_length ) );
}

/// \brief Concrete implementation that describes what this score means
string mi_score::do_description() const {
	return "(One plus the number of aligned residues"
	       + full_length_getter_ptr->description_brackets_string()
	       + ") divided by (one plus two-thirds the RMSD"
	       + rmsd.description_brackets_string()
	       + ") and by (one plus the length of the "
	       + full_length_getter_ptr->get_choice_adjective()
	       +" structure "
	       + "), all subtracted from one";
}

/// \brief TODOCUMENT
string mi_score::do_id_name() const {
	const bool length_getter_is_longer = ( full_length_getter_ptr->get_length_getter_category() == length_getter_category::LONGER );
	return "MI" + ( length_getter_is_longer ? string( "MAX" ) : string() );
}

/// \brief TODOCUMENT
str_bool_pair_vec mi_score::do_short_name_suffixes() const {
	return copy_build<str_bool_pair_vec>( join (
		length_getter_as_short_name_suffixes( *full_length_getter_ptr, { length_getter_category::LONGER, length_getter_category::SHORTER } ),
		rmsd.short_name_suffixes()
	) );
}

/// \brief Concrete implementation providing long name
string mi_score::do_long_name() const {
	return "Match Index over "
	       + full_length_getter_ptr->get_choice_adjective()
	       + rmsd.long_suffix_string();
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string mi_score::do_reference() const {
	return { "Kleywegt GJ, Jones A. Superposition. CCP4/ESF-EACBM Newsletter Protein Crystallog. 1994;31:9–14."
	         " or "
	         "Kleywegt, G. J. (1996). Use of non-crystallographic symmetry in protein structure refinement. Acta Crystallogr. D. Biol. Crystallogr.52, 842-857." };
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> mi_score::do_build_from_short_name_spec(const string &prm_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                       ) const {
//	cerr << "Should build a mi_score from string \"" << prm_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool mi_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                   ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Ctor for mi_score that allows the caller to specify the protein_only_length_getter
mi_score::mi_score(const sym_protein_only_length_getter &prm_length_choice ///< The method for choosing which of the two proteins should be used in the normalisation
                   ) : full_length_getter_ptr( prm_length_choice.sym_protein_only_clone() ) {
}

/// \brief Ctor for mi_score that allows the caller to specify the protein_only_length_getter, common_residue_selection_policy and common_atom_selection_policy
mi_score::mi_score(const sym_protein_only_length_getter  &prm_length_choice,     ///< The method for choosing which of the two proteins should be used in the normalisation
                   const common_residue_selection_policy &prm_comm_res_seln_pol, ///< The policy to use for selecting common residues
                   const common_atom_selection_policy    &prm_comm_atom_seln_pol ///< The policy to use for selecting common atoms
                   ) : rmsd              ( prm_comm_res_seln_pol, prm_comm_atom_seln_pol  ),
                       num_aligned_getter( prm_comm_res_seln_pol                          ),
                       full_length_getter_ptr( prm_length_choice.sym_protein_only_clone() ) {
}

/// \brief TODOCUMENT
const rmsd_score & mi_score::get_rmsd() const {
	return rmsd;
}

/// \brief TODOCUMENT
const num_aligned_length_getter & mi_score::get_num_aligned_getter() const {
	return num_aligned_getter;
}

/// \brief TODOCUMENT
const sym_protein_only_length_getter & mi_score::get_full_length_getter() const {
	return *full_length_getter_ptr;
}

/// \brief TODOCUMENT
///
/// \relates mi_score
bool cath::score::operator<(const mi_score &prm_mi_score_a, ///< TODOCUMENT
                            const mi_score &prm_mi_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( prm_mi_score_a, prm_mi_score_b );
	the_helper.register_comparison_field( &mi_score::get_rmsd               );
	the_helper.register_comparison_field( &mi_score::get_num_aligned_getter );
	the_helper.register_comparison_field( &mi_score::get_full_length_getter );
	return final_less_than_result( the_helper );
}
