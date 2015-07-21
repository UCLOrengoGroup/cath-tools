/// \file
/// \brief The residue_querier class definitions

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

#include "residue_querier.h"

#include <boost/numeric/conversion/cast.hpp>

#include "common/debug_numeric_cast.h"
#include "ssap/context_res.h"
#include "ssap/ssap.h"
#include "structure/protein/protein.h"

using namespace cath;
using namespace std;

using boost::numeric_cast;

/// As in the SSAP paper(s), the a and b values values are used to convert the distance into a score
/// for dynamic programming. The inherited code (this is being written in August 2013), which
/// appears to use the square of the distance between residues rather than the distance as indicated
/// in the SSAP paper. Otherwise, the formula for the score s the same: \f$ s = \frac{a}{b + d^2} \f$

/// \brief The value a used in the SSAP paper (for residues)
/////
///// Note that this scaled by INTEGER_SCALING^2 ( = 10 * 10 = 100), which matches
///// the result of the distance being scaled by INTEGER_SCALING
//const size_t residue_querier::RESIDUE_A_VALUE          = 500 * entry_querier::INTEGER_SCALING * entry_querier::INTEGER_SCALING;
const float_score_type residue_querier::RESIDUE_A_VALUE          = 500.0;

/// \brief The value b used in the SSAP paper (for residues)
/////
///// Note that this scaled by INTEGER_SCALING^2 ( = 10 * 10 = 100), which matches
///// the result of the distance being scaled by INTEGER_SCALING
//const size_t residue_querier::RESIDUE_B_VALUE          =  10 * entry_querier::INTEGER_SCALING * entry_querier::INTEGER_SCALING;
const float_score_type residue_querier::RESIDUE_B_VALUE          =  10.0;

/// \brief The minimum score that the algorithm will consider for residues
///        (often referred to as c)
///
/// Note that this does not need scaling
const float_score_type residue_querier::RESIDUE_MIN_SCORE_CUTOFF =  10.0;

/// \brief The maximum squared-distance value that the algorithm will consider (for residues)
///
/// This is calculated as \f$ \frac{a}{c} - b\f$,
/// which comes from solving \f$ c = \frac{a}{ d^2 + b} \f$ for \f$ s^2 \f$ )
///
/// NOTE: This has an implied scaling ( * INTEGER_SCALING * INTEGER_SCALING) built into it.
///       This line will not need changing as long as A, B and the distances have their scaling
///       removed simultaneously.
///
/// With current (unscaled) values of a=500, b=10 and c=10, this gives a maximum (unscaled) value
/// of 40 for d^2, which implies a maximum value for d of sqrt(40) ~= 6.3A
const float_score_type residue_querier::RESIDUE_MAX_DIST_SQ_CUTOFF(
	residue_querier::RESIDUE_A_VALUE / residue_querier::RESIDUE_MIN_SCORE_CUTOFF - residue_querier::RESIDUE_B_VALUE
);

/// \brief TODOCUMENT
size_t residue_querier::do_get_length(const protein &arg_protein ///< TODOCUMENT
                                      ) const {
	return arg_protein.get_length();
}

/// \brief TODOCUMENT
double residue_querier::do_get_gap_penalty_ratio() const {
	return 1.0;
}

/// \brief TODOCUMENT
size_t residue_querier::do_num_excluded_on_either_size() const {
	return 5;
}

/// \brief TODOCUMENT
float_score_type residue_querier::do_optimum_single_score() const {
	return RESIDUE_A_VALUE / RESIDUE_B_VALUE;
}

/// \brief TODOCUMENT
string residue_querier::do_get_entry_name() const {
	return "residue";
}

/// \brief TODOCUMENT
///
/// TEMPORARILY DEBUGGING
score_type residue_querier::do_distance_score__offset_1(const protein &arg_protein_a,                   ///< TODOCUMENT
                                                        const protein &arg_protein_b,                   ///< TODOCUMENT
                                                        const size_t  &arg_a_view_from_index__offset_1, ///< TODOCUMENT
                                                        const size_t  &arg_b_view_from_index__offset_1, ///< TODOCUMENT
                                                        const size_t  &arg_a_dest_to_index__offset_1,   ///< TODOCUMENT
                                                        const size_t  &arg_b_dest_to_index__offset_1    ///< TODOCUMENT
                                                        ) const {
	const residue &residue_a_view_from = get_residue_ref_of_index__offset_1( arg_protein_a, arg_a_view_from_index__offset_1 );
	const residue &residue_b_view_from = get_residue_ref_of_index__offset_1( arg_protein_b, arg_b_view_from_index__offset_1 );
	const residue &residue_a_dest_to   = get_residue_ref_of_index__offset_1( arg_protein_a, arg_a_dest_to_index__offset_1   );
	const residue &residue_b_dest_to   = get_residue_ref_of_index__offset_1( arg_protein_b, arg_b_dest_to_index__offset_1   );
	return debug_numeric_cast<score_type>(
		context_res<true>(
			residue_a_view_from, residue_b_view_from,
			residue_a_dest_to,   residue_b_dest_to
		)
	);
}

/// \brief TODOCUMENT
bool residue_querier::do_are_comparable__offset_1(const protein &/*arg_protein_a*/,                   ///< TODOCUMENT
                                                  const protein &/*arg_protein_b*/,                   ///< TODOCUMENT
                                                  const size_t  &/*arg_a_view_from_index__offset_1*/, ///< TODOCUMENT
                                                  const size_t  &/*arg_b_view_from_index__offset_1*/, ///< TODOCUMENT
                                                  const size_t  &/*arg_a_dest_to_index__offset_1*/,   ///< TODOCUMENT
                                                  const size_t  &/*arg_b_dest_to_index__offset_1*/    ///< TODOCUMENT
                                                  ) const {
	return true;
}

/// \brief TODOCUMENT
bool residue_querier::do_are_similar__offset_1(const protein &arg_protein_a,         ///< TODOCUMENT
                                               const protein &arg_protein_b,         ///< TODOCUMENT
                                               const size_t  &arg_index_a__offset_1, ///< TODOCUMENT
                                               const size_t  &arg_index_b__offset_1  ///< TODOCUMENT
                                               ) const {
	const residue &residue_a = get_residue_ref_of_index__offset_1( arg_protein_a, arg_index_a__offset_1 );
	const residue &residue_b = get_residue_ref_of_index__offset_1( arg_protein_b, arg_index_b__offset_1 );
	return residues_have_similar_area_angle_props(residue_a, residue_b);
}

/// \brief TODOCUMENT
bool residue_querier::do_temp_hacky_is_residue() const {
	return true;
}
