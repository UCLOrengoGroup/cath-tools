/// \file
/// \brief The sequence_similarity_score class definitions

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

#include "sequence_similarity_score.hpp"

//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/archive/xml_oarchive.hpp>
#include <boost/logic/tribool.hpp>
#include <boost/range/join.hpp>
//#include <boost/serialization/export.hpp>

#include "alignment/alignment.hpp"
#include "alignment/pair_alignment.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/less_than_helper.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "score/length_getter/length_of_longer_getter.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/protein/amino_acid.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/structure_type_aliases.hpp"
#include "superposition/superposition.hpp"

using namespace boost::logic;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::score;
using namespace cath::sup;
using namespace std;

using boost::numeric_cast;
using boost::range::join;

//BOOST_CLASS_EXPORT(sequence_similarity_score)

/// \brief A standard do_clone method.
unique_ptr<aligned_pair_score> sequence_similarity_score::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Concrete implementation that records that a lower RMSD score is better
tribool sequence_similarity_score::do_higher_is_better() const {
	return true;
}

/// \brief Concrete implementation for calculating the RMSD of an alignment
score_value sequence_similarity_score::do_calculate(const alignment &arg_alignment, ///< The pair alignment to be scored
                                                    const protein   &arg_protein_a, ///< The protein associated with the first  half of the alignment
                                                    const protein   &arg_protein_b  ///< The protein associated with the second half of the alignment
                                                    ) const {
	const size_t length = arg_alignment.length();

	score_type score( 0 );
	for (const size_t &index : indices( length ) ) {
		if ( has_both_positions_of_index( arg_alignment, index ) ) {
			const aln_posn_type a_posn  = get_a_position_of_index( arg_alignment, index  );
			const aln_posn_type b_posn  = get_b_position_of_index( arg_alignment, index  );
			const amino_acid    amino_a = get_amino_acid_of_index( arg_protein_a, a_posn );
			const amino_acid    amino_b = get_amino_acid_of_index( arg_protein_b, b_posn );
			score += scores.get_score( amino_a, amino_b );
		}
	}

	const size_t      normalisation_length = length_getter_ptr->get_length( arg_alignment, arg_protein_a, arg_protein_b );
	const score_value normalisation_score  = numeric_cast<score_value>( normalisation_length * numeric_cast<size_t>( scores.get_highest_score() ) );
	return 100.0 * numeric_cast<score_value>( score ) / normalisation_score;
}

/// \brief Concrete implementation that describes what this score means
string sequence_similarity_score::do_description() const {
	return "The sequence similarity between the two aligned sequences, as calculated using a " + scores.get_name() + " matrix"
	       + length_getter_ptr->description_brackets_string();
}

/// \brief TODOCUMENT
string sequence_similarity_score::do_id_name() const {
	return "sequence_id";
}

/// \brief TODOCUMENT
str_bool_pair_vec sequence_similarity_score::do_short_name_suffixes() const {
	const auto subs_matrix_suffixes = { make_pair( scores.get_name(), true ) };
	return copy_build<str_bool_pair_vec>( join (
		subs_matrix_suffixes,
		length_getter_as_short_name_suffixes( *length_getter_ptr, { length_getter_category::LONGER } )
	) );
}

/// \brief Concrete implementation providing long name
string sequence_similarity_score::do_long_name() const {
	return "Sequence identity (using " + scores.get_name() + ", normalised over " + length_getter_ptr->long_name() + ")";
}

///// \brief Build an aligned_pair_score of this concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> sequence_similarity_score::do_build_from_short_name_spec(const string &arg_short_name_spec ///< The short_name_spec that defines any properties that the resulting aligned_pair_score should have
//                                                                                        ) const {
//	cerr << "Should build a sequence_similarity_score from string \"" << arg_short_name_spec << "\"" << endl;
//	return clone();
//}

/// \brief TODOCUMENT
bool sequence_similarity_score::do_less_than_with_same_dynamic_type(const aligned_pair_score &arg_aligned_pair_score ///< TODOCUMENT
                                                                    ) const {
	const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( arg_aligned_pair_score );
	return ( *this < casted_aligned_pair_score );
}

/// \brief Default ctor for sequence_similarity_score (the substitution matrix defaults to the identity matrix)
sequence_similarity_score::sequence_similarity_score(substitution_matrix arg_scores ///< The substitution matrix to use for scoring
                                                     ) : scores { std::move( arg_scores ) } {
}

/// \brief Ctor for sequence_similarity_score that allows the caller to specify the protein_only_length_getter
///        (the substitution matrix defaults to the identity matrix)
sequence_similarity_score::sequence_similarity_score(const length_getter &arg_length_getter, ///< The method for choosing the length used for normalisation
                                                     substitution_matrix  arg_scores         ///< The substitution matrix to use for scoring
                                                     ) : scores            { std::move( arg_scores )   },
                                                         length_getter_ptr { arg_length_getter.clone() } {
}

/// \brief TODOCUMENT
///
/// \relates sequence_similarity_score
const substitution_matrix & sequence_similarity_score::get_substitution_matrix() const {
	return scores;
}

/// \brief TODOCUMENT
///
/// \relates sequence_similarity_score
const length_getter & sequence_similarity_score::get_length_getter() const {
	return *length_getter_ptr;
}

/// \brief TODOCUMENT
///
/// \relates sequence_similarity_score
bool cath::score::operator<(const sequence_similarity_score &arg_sequence_similarity_score_a, ///< TODOCUMENT
                            const sequence_similarity_score &arg_sequence_similarity_score_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( arg_sequence_similarity_score_a, arg_sequence_similarity_score_b );
	the_helper.register_comparison_field( &sequence_similarity_score::get_substitution_matrix );
	the_helper.register_comparison_field( &sequence_similarity_score::get_length_getter       );
	return final_less_than_result( the_helper );
}
