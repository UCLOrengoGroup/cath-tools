/// \file
/// \brief The sec_struc_querier class definitions

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

#include "sec_struc_querier.h"

#include <boost/numeric/conversion/cast.hpp>

#include "ssap/ssap.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;

/// \brief The value a used in the SSAP paper (for secondary structures)
///
/// Note that this scaled by INTEGER_SCALING^2 ( = 10 * 10 = 100), which matches
/// the result of the distance being scaled by INTEGER_SCALING
const size_t sec_struc_querier::SEC_STRUC_A_VALUE          = 800 * entry_querier::INTEGER_SCALING * entry_querier::INTEGER_SCALING;

/// \brief The value b used in the SSAP paper (for secondary structures)
///
/// Note that this scaled by INTEGER_SCALING^2 ( = 10 * 10 = 100), which matches
/// the result of the distance being scaled by INTEGER_SCALING
const size_t sec_struc_querier::SEC_STRUC_B_VALUE          =   6 * entry_querier::INTEGER_SCALING * entry_querier::INTEGER_SCALING;

/// \brief The minimum score that the algorithm will consider for secondary structures
///        (often referred to as c)
///
/// Note that this does not need scaling
const size_t sec_struc_querier::SEC_STRUC_MIN_SCORE_CUTOFF =  10;

/// \brief The maximum squared-distance value that the algorithm will consider (for secondary structures)
///
/// This is calculated as \f$ \frac{a}{c} - b\f$,
/// which comes from solving \f$ c = \frac{a}{ d^2 + b} \f$ for \f$ s^2 \f$ )
///
/// NOTE: This has an implied scaling ( * INTEGER_SCALING * INTEGER_SCALING) built into it.
///       This line will not need changing as long as A, B and the distances have their scaling
///       removed simultaneously.
///
/// With current (unscaled) values of a=800, b=6 and c=10, this gives a maximum (unscaled) value
/// of 74 for d^2, which implies a maximum value for d of sqrt(74) ~= 8.6A
const size_t sec_struc_querier::SEC_STRUC_MAX_DIST_SQ_CUTOFF(
	sec_struc_querier::SEC_STRUC_A_VALUE / sec_struc_querier::SEC_STRUC_MIN_SCORE_CUTOFF - sec_struc_querier::SEC_STRUC_B_VALUE
);

/// \brief TODOCUMENT
size_t sec_struc_querier::do_get_length(const protein &arg_protein ///< TODOCUMENT
                                        ) const {
	return arg_protein.get_num_sec_strucs();
}

/// \brief TODOCUMENT
double sec_struc_querier::do_get_gap_penalty_ratio() const {
	return 0.0375;
}

/// \brief TODOCUMENT
size_t sec_struc_querier::do_num_excluded_on_either_size() const {
	return 0;
}

/// \brief TODOCUMENT
double sec_struc_querier::do_optimum_single_score() const {
	return numeric_cast<double>(SEC_STRUC_A_VALUE) / numeric_cast<double>(SEC_STRUC_B_VALUE);
}

/// \brief TODOCUMENT
string sec_struc_querier::do_get_entry_name() const {
	return "secondary structure";
}

/// \brief TODOCUMENT
score_type sec_struc_querier::do_distance_score__offset_1(const protein &arg_protein_a,                   ///< TODOCUMENT
                                                          const protein &arg_protein_b,                   ///< TODOCUMENT
                                                          const size_t  &arg_a_view_from_index__offset_1, ///< TODOCUMENT
                                                          const size_t  &arg_b_view_from_index__offset_1, ///< TODOCUMENT
                                                          const size_t  &arg_a_dest_to_index__offset_1,   ///< TODOCUMENT
                                                          const size_t  &arg_b_dest_to_index__offset_1    ///< TODOCUMENT
                                                          ) const {
	// Sanity check the inputs
	const size_t num_sec_strucs_in_a = arg_protein_a.get_num_sec_strucs();
	if (arg_a_view_from_index__offset_1 < 1 || arg_a_view_from_index__offset_1 > num_sec_strucs_in_a ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("arg_a_view_from_index__offset_1 is out of range"));
	}
	if (arg_a_dest_to_index__offset_1   < 1 || arg_a_dest_to_index__offset_1   > num_sec_strucs_in_a ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("arg_a_dest_to_index__offset_1   is out of range"));
	}
	const size_t num_sec_strucs_in_b = arg_protein_b.get_num_sec_strucs();
	if (arg_b_view_from_index__offset_1 < 1 || arg_b_view_from_index__offset_1 > num_sec_strucs_in_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("arg_b_view_from_index__offset_1 is out of range"));
	}
	if (arg_b_dest_to_index__offset_1   < 1 || arg_b_dest_to_index__offset_1   > num_sec_strucs_in_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("arg_b_dest_to_index__offset_1   is out of range"));
	}

	// Pass through to context_sec (whilst switching the indices to use offset 0)
	return context_sec(
		arg_protein_a,                       arg_protein_b,
		arg_a_view_from_index__offset_1 - 1, arg_b_view_from_index__offset_1 - 1,
		arg_a_dest_to_index__offset_1   - 1, arg_b_dest_to_index__offset_1   - 1
	);
}

/// \brief TODOCUMENT
bool sec_struc_querier::do_are_comparable__offset_1(const protein &arg_protein_a,                   ///< TODOCUMENT
                                                    const protein &arg_protein_b,                   ///< TODOCUMENT
                                                    const size_t  &arg_a_view_from_index__offset_1, ///< TODOCUMENT
                                                    const size_t  &arg_b_view_from_index__offset_1, ///< TODOCUMENT
                                                    const size_t  &arg_a_dest_to_index__offset_1,   ///< TODOCUMENT
                                                    const size_t  &arg_b_dest_to_index__offset_1    ///< TODOCUMENT
                                                    ) const {
	check_offset_1( arg_a_view_from_index__offset_1 );
	check_offset_1( arg_b_view_from_index__offset_1 );
	check_offset_1( arg_a_dest_to_index__offset_1   );
	check_offset_1( arg_b_dest_to_index__offset_1   );

	const sec_struc &sec_struc_i_a = arg_protein_a.get_sec_struc_ref_of_index( arg_a_view_from_index__offset_1 - 1 );
	const sec_struc &sec_struc_i_b = arg_protein_b.get_sec_struc_ref_of_index( arg_b_view_from_index__offset_1 - 1 );
	const sec_struc &sec_struc_j_a = arg_protein_a.get_sec_struc_ref_of_index( arg_a_dest_to_index__offset_1   - 1 );
	const sec_struc &sec_struc_j_b = arg_protein_b.get_sec_struc_ref_of_index( arg_b_dest_to_index__offset_1   - 1 );

	const bool i_sec_strucs_match  = ( sec_struc_i_a.get_type() != sec_struc_i_b.get_type() );
	const bool j_sec_strucs_match  = ( sec_struc_j_a.get_type() != sec_struc_j_b.get_type() );

	return ( i_sec_strucs_match && j_sec_strucs_match );
}

/// \brief TODOCUMENT
bool sec_struc_querier::do_are_similar__offset_1(const protein &arg_protein_a,         ///< TODOCUMENT
                                                 const protein &arg_protein_b,         ///< TODOCUMENT
                                                 const size_t  &arg_index_a__offset_1, ///< TODOCUMENT
                                                 const size_t  &arg_index_b__offset_1  ///< TODOCUMENT
                                                 ) const {
	check_offset_1( arg_index_a__offset_1 );
	check_offset_1( arg_index_b__offset_1 );

	const sec_struc &sec_struc_a = arg_protein_a.get_sec_struc_ref_of_index( arg_index_a__offset_1 - 1 );
	const sec_struc &sec_struc_b = arg_protein_b.get_sec_struc_ref_of_index( arg_index_b__offset_1 - 1 );

	return ( sec_struc_a.get_type() == sec_struc_b.get_type() );
}

/// \brief TODOCUMENT
bool sec_struc_querier::do_temp_hacky_is_residue() const {
	return false;
}

