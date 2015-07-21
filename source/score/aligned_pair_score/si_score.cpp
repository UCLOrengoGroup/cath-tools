/// \file
/// \brief The si_score class definitions

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

#include "si_score.h"

//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/archive/xml_oarchive.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/range/join.hpp>
//#include <boost/serialization/export.hpp>

#include "alignment/alignment.h"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.h"
#include "common/algorithm/copy_build.h"
#include "common/clone/make_uptr_clone.h"
#include "common/less_than_helper.h"
#include "score/length_getter/length_of_longer_getter.h"
#include "score/length_getter/length_of_shorter_getter.h"

#include <iostream> // ***** TEMPORARY *****

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::range::join;

//BOOST_CLASS_EXPORT(si_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> si_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower SAS is generally better
tribool si_score::do_higher_is_better() const {
	return false;
}

/// \brief Concrete implementation for calculating the SAS of an alignment
score_value si_score::do_calculate(const alignment &arg_alignment, ///< The pair alignment to be scored
                                   const protein   &arg_protein_a, ///< The protein associated with the first  half of the alignment
                                   const protein   &arg_protein_b  ///< The protein associated with the second half of the alignment
                                   ) const {
	return rmsd.calculate( arg_alignment, arg_protein_a, arg_protein_b )
	       * get_length_score( *full_length_getter_ptr,                arg_protein_a, arg_protein_b )
	       / get_length_score(  num_aligned_getter,     arg_alignment, arg_protein_a, arg_protein_b );
}

/// \brief Concrete implementation that describes what this score means
string si_score::do_description() const {
	return "The RMSD"
	       + rmsd.description_brackets_string()
	       + " divided by the fraction of the "
	       + full_length_getter_ptr->get_choice_adjective()
	       +" structure's residues that have been aligned"
	       + num_aligned_getter.description_brackets_string();
}

/// \brief TODOCUMENT
string si_score::do_id_name() const {
	const bool length_getter_is_longer = ( full_length_getter_ptr->get_length_getter_category() == length_getter_category::LONGER );
	return "SI" + ( length_getter_is_longer ? string( "MAX" ) : string() );
}

/// \brief TODOCUMENT
str_bool_pair_vec si_score::do_short_name_suffixes() const {
	return copy_build<str_bool_pair_vec>( join (
		length_getter_as_short_name_suffixes( *full_length_getter_ptr, { length_getter_category::LONGER, length_getter_category::SHORTER } ),
		rmsd.short_name_suffixes()
	) );
}

/// \brief Concrete implementation providing long name
string si_score::do_long_name() const {
	return "Similarity Index over "
	       + full_length_getter_ptr->get_choice_adjective()
	       + rmsd.long_suffix_string();
}

/// \brief Concrete implementation for providing a reference to a publication describing this score
string si_score::do_reference() const {
	return { "Kleywegt GJ, Jones A. Superposition. CCP4/ESF-EACBM Newsletter Protein Crystallog. 1994;31:9–14."
	         " or "
	         "Kleywegt, G. J. (1996). Use of non-crystallographic symmetry in protein structure refinement. Acta Crystallogr. D. Biol. Crystallogr.52, 842-857." };
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> si_score::do_build_from_short_name_spec(const string &arg_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                       ) const {
//	cerr << "Should build a si_score from string \"" << arg_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool si_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &arg_aligned_pair_score ///< TODOCUMENT
                                                   ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( arg_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Ctor for si_score that allows the caller to specify the protein_only_length_getter
si_score::si_score(const sym_protein_only_length_getter &arg_length_choice ///< The method for choosing which of the two proteins should be used in the normalisation
                   ) : full_length_getter_ptr( arg_length_choice.sym_protein_only_clone() ) {
}

/// \brief Ctor for si_score that allows the caller to specify the protein_only_length_getter, common_residue_selection_policy and common_atom_selection_policy
si_score::si_score(const sym_protein_only_length_getter  &arg_length_choice,     ///< The method for choosing which of the two proteins should be used in the normalisation
                   const common_residue_selection_policy &arg_comm_res_seln_pol, ///< The policy to use for selecting common residues
                   const common_atom_selection_policy    &arg_comm_atom_seln_pol ///< The policy to use for selecting common atoms
                   ) : rmsd                  ( arg_comm_res_seln_pol, arg_comm_atom_seln_pol ),
                       num_aligned_getter    ( arg_comm_res_seln_pol                         ),
                       full_length_getter_ptr( arg_length_choice.sym_protein_only_clone()    ) {
}

/// \brief TODOCUMENT
const rmsd_score & si_score::get_rmsd() const {
	return rmsd;
}

/// \brief TODOCUMENT
const num_aligned_length_getter & si_score::get_num_aligned_getter() const {
	return num_aligned_getter;
}

/// \brief TODOCUMENT
const sym_protein_only_length_getter & si_score::get_full_length_getter() const {
	return *full_length_getter_ptr;
}

/// \brief TODOCUMENT
///
/// \relates si_score
si_score cath::score::make_simax_score(const common_residue_selection_policy &arg_common_residue_selection_policy, ///< TODOCUMENT
                                       const common_atom_selection_policy    &arg_common_atom_selection_policy     ///< TODOCUMENT
                                       ) {
	return si_score(
		length_of_longer_getter(),
		arg_common_residue_selection_policy,
		arg_common_atom_selection_policy
	);
}

/// \brief TODOCUMENT
///
/// \relates si_score
bool cath::score::operator<(const si_score &arg_si_score_a, ///< TODOCUMENT
                            const si_score &arg_si_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( arg_si_score_a, arg_si_score_b );
	the_helper.register_comparison_field( &si_score::get_rmsd               );
	the_helper.register_comparison_field( &si_score::get_num_aligned_getter );
	the_helper.register_comparison_field( &si_score::get_full_length_getter );
	return final_less_than_result( the_helper );
}
