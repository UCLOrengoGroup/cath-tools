/// \file
/// \brief The substitution_matrix class definitions

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

#include "substitution_matrix.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm.hpp>

#include "common/algorithm/sort_uniq_copy.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/invert_permutation.hpp"
#include "common/less_than_helper.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/not_implemented_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "score/aligned_pair_score/substitution_matrix/blosum62_substitution_matrix.hpp"
#include "score/aligned_pair_score/substitution_matrix/identity_substitution_matrix.hpp"
#include "score/aligned_pair_score/substitution_matrix/match_substitution_matrix.hpp"
#include "structure/protein/amino_acid.hpp"
#include "structure/structure_type_aliases.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::numeric_cast;
using boost::range::max_element;

/// \brief TODOCUMENT
const amino_acid_vec & substitution_matrix::get_amino_acids() const {
	return amino_acids;
}

/// \brief TODOCUMENT
const score_vec_vec & substitution_matrix::get_scores() const {
	return scores;
}

/// \brief TODOCUMENT
const score_type & substitution_matrix::get_score_for_one_unknown_aa() const {
	return score_for_one_unknown_aa;
}

/// \brief TODOCUMENT
const score_type & substitution_matrix::get_score_for_two_unknown_aas() const {
	return score_for_two_unknown_aas;
}

/// \brief Construct the permutation represented by changing from the original amino acid ordering to the new amino acid ordering
///
/// \todo Tidied up this code
size_vec substitution_matrix::order_permutation(const amino_acid_vec &arg_orig_aa_ordering, ///< The original amino acid ordering
                                                const amino_acid_vec &arg_new_aa_ordering   ///< The new amino acid ordering
                                                ) {
	// Grab the number of amino acids and check it's consistent
	const size_t num_amino_acids = arg_orig_aa_ordering.size();
	if ( num_amino_acids != arg_new_aa_ordering.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Mismatch in the number of amino acids"));
	}

	// Construct a permutation and reserve the memory
	size_vec permutation;
	permutation.reserve( num_amino_acids );

	// For each of the amino acids in the original ordering...
	for (const amino_acid &orig_amino_acid : arg_orig_aa_ordering) {
		// Find the index of the amino acid in the new ordering
		const size_t index_in_new = numeric_cast<size_t>( distance(
			common::cbegin( arg_new_aa_ordering ),
			find(
				arg_new_aa_ordering,
				orig_amino_acid
			)
		) );
		// ...check the index is valid...
		if ( index_in_new >= num_amino_acids ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid not found in new amino acids"));
		}
		// ...and then add the index to the back of the permutation
		permutation.push_back( index_in_new );
	}
	// Return the resulting permutation
	return permutation;
}

/// \brief Reorder the scores corresponding to the first amino acid ordering to the second amino acid ordering
score_vec_vec substitution_matrix::reorder_scores(const amino_acid_vec &arg_orig_aa_ordering, ///< The original amino acid ordering
                                                  const amino_acid_vec &arg_new_aa_ordering,  ///< The required, new amino acid ordering
                                                  const score_vec_vec  &arg_scores            ///< Scores corresponding to the original amino acid ordering
                                                  ) {
	// Get the permutation required to get from the new ordering back to the old one
	const size_vec back_permutation = invert_permutation(
		order_permutation( arg_orig_aa_ordering, arg_new_aa_ordering )
	);

	// Construct a matrix to hold the new scores
	const size_t num_amino_acids = back_permutation.size();
	score_vec_vec scores( num_amino_acids, score_vec( num_amino_acids, 0 ) );

	// Loop over the new matrix, populating its scores from the old matrix
	for (size_t new_aa_index_a = 0; new_aa_index_a < num_amino_acids; ++new_aa_index_a) {
		const size_t orig_aa_index_a = back_permutation[ new_aa_index_a ];
		for (size_t new_aa_index_b = 0; new_aa_index_b < num_amino_acids; ++new_aa_index_b) {
			const size_t orig_aa_index_b = back_permutation[ new_aa_index_b ];
			scores[ new_aa_index_a ][ new_aa_index_b ] = arg_scores[ orig_aa_index_a ][ orig_aa_index_b ];
		}
	}

	// Return the results
	return scores;
}

/// \brief Return whether one row/column has a lower highest score than the next
bool substitution_matrix::has_lower_highest_score(const score_vec &arg_scores_a, ///< The first row/column of scores to compare
                                                  const score_vec &arg_scores_b  ///< The second row/column of scores to compare
                                                  ) {
	return ( highest_score_of_scores( arg_scores_a ) < highest_score_of_scores( arg_scores_b ) );
}

/// \brief Return the highest score to be found in the specified row/column
score_type substitution_matrix::highest_score_of_scores(const score_vec &arg_scores ///< The row/column of scores to examine
                                                        ) {
	return *max_element( arg_scores );
}

/// \brief Check that the matrix is symmetric
///
/// \pre The matrix must be symmetric else an out_of_range_exception will be thrown
void substitution_matrix::check_is_symmetric() const {
	const size_t num_amino_acids = amino_acids.size();
	for (size_t index_a = 0; index_a < num_amino_acids; ++index_a) {
		for (size_t index_b = 0; index_b < num_amino_acids; ++index_b) {
			if ( scores[ index_a ][ index_b ] != scores[ index_b ][ index_a ] ) {
				BOOST_THROW_EXCEPTION(out_of_range_exception("The substitution_matrix data isn't symmetric"));
			}
		}
	}
}

/// \brief Ctor for substitution_matrix
///
/// arg_amino_acids needn't necessarily be sorted
substitution_matrix::substitution_matrix(const amino_acid_vec &arg_amino_acids,               ///< The list of amino acids
                                         const score_vec_vec  &arg_scores,                    ///< The all-against-all scores corresponding to the list of amino acids
                                         const score_type     &arg_score_for_one_unknown_aa,  ///< The score that's returned if one of the query amino acids is unrecognised
                                         const score_type     &arg_score_for_two_unknown_aas, ///< The score that's returned if both query amino acids are unrecognised
                                         const string         &arg_name                       ///< A name for the substitution matrix
                                         ) : amino_acids              ( sort_copy( arg_amino_acids )  ),
                                             scores                   ( reorder_scores( arg_amino_acids, amino_acids, arg_scores ) ),
                                             score_for_one_unknown_aa ( arg_score_for_one_unknown_aa  ),
                                             score_for_two_unknown_aas( arg_score_for_two_unknown_aas ),
                                             name                     ( arg_name                      ) {
	// Check that the generated matrix is symmetric
	check_is_symmetric();

	if ( scores.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot create empty substitution matrix"));
	}
}

/// \brief Getter for the name
const string & substitution_matrix::get_name() const {
	return name;
}

/// \brief The highest score of any of the elements in the matrix
score_type substitution_matrix::get_highest_score() const {
	return highest_score_of_scores(
		*max_element(
			scores,
			&has_lower_highest_score
		)
	);
}

/// \brief Get the score corresponding to the two specified amino acids
///
/// \todo Tidied up this code
score_type substitution_matrix::get_score(const amino_acid &arg_amino_acid_a, ///< The first amino acid of the query
                                          const amino_acid &arg_amino_acid_b  ///< The second amino acid of the query
                                          ) const {
	// Grab the index of the first amino acid in amino_acids
	const size_t index_a = numeric_cast<size_t>( distance(
		common::cbegin( amino_acids ),
		lower_bound( amino_acids, arg_amino_acid_a )
	) );

	// Grab the index of the second amino acid in amino_acids
	const size_t index_b = numeric_cast<size_t>( distance(
		common::cbegin( amino_acids ),
		lower_bound( amino_acids, arg_amino_acid_b )
	) );

	// Bounds check both results
	const size_t num_amino_acids = amino_acids.size();
	if ( index_a > num_amino_acids || index_b > num_amino_acids ) {
		const bool both_unknown = ( index_a > num_amino_acids && index_b > num_amino_acids );
		return both_unknown ? score_for_two_unknown_aas
		                    : score_for_one_unknown_aa;
	}

	// Return the relevant score
	return scores[ index_a ][ index_b ];
}

/// \brief TODOCUMENT
///
/// \relates substitution_matrix
bool cath::score::operator<(const substitution_matrix &arg_substitution_matrix_a, ///< TODOCUMENT
                            const substitution_matrix &arg_substitution_matrix_b  ///< TODOCUMENT
                            ) {
	auto the_helper = make_less_than_helper( arg_substitution_matrix_a, arg_substitution_matrix_b );

	the_helper.register_comparison_field( &substitution_matrix::get_name                      );
	the_helper.register_comparison_field( &substitution_matrix::get_amino_acids               );
	the_helper.register_comparison_field( &substitution_matrix::get_scores                    );
	the_helper.register_comparison_field( &substitution_matrix::get_score_for_one_unknown_aa  );
	the_helper.register_comparison_field( &substitution_matrix::get_score_for_two_unknown_aas );

	return final_less_than_result( the_helper );
}

/// \brief Make a standard list of amino acids
///
/// \relates substitution_matrix
amino_acid_vec cath::score::get_standard_substitution_amino_acids() {
	return make_amino_acids_of_chars( {
		'A', 'B', 'C', 'D', 'E',
		'F', 'G', 'H', 'I', 'K',
		'L', 'M', 'N', 'P', 'Q',
		'R', 'S', 'T', 'V', 'W',
		'X', 'Y', 'Z'
	} );
}

/// \brief Get a list of all the
///
/// \relates substitution_matrix
substitution_matrix_vec cath::score::get_all_substitution_matrices() {
	return {
		make_subs_matrix_blosum62(),
		make_subs_matrix_identity(),
		make_subs_matrix_match()
	};
}
