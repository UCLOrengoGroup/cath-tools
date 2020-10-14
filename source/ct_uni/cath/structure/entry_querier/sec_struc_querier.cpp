/// \file
/// \brief The sec_struc_querier class definitions

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

#include "sec_struc_querier.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "cath/ssap/ssap.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;

constexpr size_t sec_struc_querier::SEC_STRUC_A_VALUE;
constexpr size_t sec_struc_querier::SEC_STRUC_B_VALUE;
constexpr size_t sec_struc_querier::SEC_STRUC_MIN_SCORE_CUTOFF;
constexpr size_t sec_struc_querier::SEC_STRUC_MAX_DIST_SQ_CUTOFF;

/// \brief TODOCUMENT
size_t sec_struc_querier::do_get_length(const protein &prm_protein ///< TODOCUMENT
                                        ) const {
	return prm_protein.get_num_sec_strucs();
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
score_type sec_struc_querier::do_distance_score__offset_1(const protein &prm_protein_a,                   ///< TODOCUMENT
                                                          const protein &prm_protein_b,                   ///< TODOCUMENT
                                                          const size_t  &prm_a_view_from_index__offset_1, ///< TODOCUMENT
                                                          const size_t  &prm_b_view_from_index__offset_1, ///< TODOCUMENT
                                                          const size_t  &prm_a_dest_to_index__offset_1,   ///< TODOCUMENT
                                                          const size_t  &prm_b_dest_to_index__offset_1    ///< TODOCUMENT
                                                          ) const {
	// Sanity check the inputs
	const size_t num_sec_strucs_in_a = prm_protein_a.get_num_sec_strucs();
	if (prm_a_view_from_index__offset_1 < 1 || prm_a_view_from_index__offset_1 > num_sec_strucs_in_a ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("prm_a_view_from_index__offset_1 is out of range"));
	}
	if (prm_a_dest_to_index__offset_1   < 1 || prm_a_dest_to_index__offset_1   > num_sec_strucs_in_a ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("prm_a_dest_to_index__offset_1   is out of range"));
	}
	const size_t num_sec_strucs_in_b = prm_protein_b.get_num_sec_strucs();
	if (prm_b_view_from_index__offset_1 < 1 || prm_b_view_from_index__offset_1 > num_sec_strucs_in_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("prm_b_view_from_index__offset_1 is out of range"));
	}
	if (prm_b_dest_to_index__offset_1   < 1 || prm_b_dest_to_index__offset_1   > num_sec_strucs_in_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("prm_b_dest_to_index__offset_1   is out of range"));
	}

	// Pass through to context_sec (whilst switching the indices to use offset 0)
	return context_sec(
		prm_protein_a,                       prm_protein_b,
		prm_a_view_from_index__offset_1 - 1, prm_b_view_from_index__offset_1 - 1,
		prm_a_dest_to_index__offset_1   - 1, prm_b_dest_to_index__offset_1   - 1
	);
}

/// \brief TODOCUMENT
bool sec_struc_querier::do_are_comparable__offset_1(const protein &prm_protein_a,                   ///< TODOCUMENT
                                                    const protein &prm_protein_b,                   ///< TODOCUMENT
                                                    const size_t  &prm_a_view_from_index__offset_1, ///< TODOCUMENT
                                                    const size_t  &prm_b_view_from_index__offset_1, ///< TODOCUMENT
                                                    const size_t  &prm_a_dest_to_index__offset_1,   ///< TODOCUMENT
                                                    const size_t  &prm_b_dest_to_index__offset_1    ///< TODOCUMENT
                                                    ) const {
	check_offset_1( prm_a_view_from_index__offset_1 );
	check_offset_1( prm_b_view_from_index__offset_1 );
	check_offset_1( prm_a_dest_to_index__offset_1   );
	check_offset_1( prm_b_dest_to_index__offset_1   );

	const sec_struc &sec_struc_i_a = prm_protein_a.get_sec_struc_ref_of_index( prm_a_view_from_index__offset_1 - 1 );
	const sec_struc &sec_struc_i_b = prm_protein_b.get_sec_struc_ref_of_index( prm_b_view_from_index__offset_1 - 1 );
	const sec_struc &sec_struc_j_a = prm_protein_a.get_sec_struc_ref_of_index( prm_a_dest_to_index__offset_1   - 1 );
	const sec_struc &sec_struc_j_b = prm_protein_b.get_sec_struc_ref_of_index( prm_b_dest_to_index__offset_1   - 1 );

	const bool i_sec_strucs_match  = ( sec_struc_i_a.get_type() != sec_struc_i_b.get_type() );
	const bool j_sec_strucs_match  = ( sec_struc_j_a.get_type() != sec_struc_j_b.get_type() );

	return ( i_sec_strucs_match && j_sec_strucs_match );
}

/// \brief TODOCUMENT
bool sec_struc_querier::do_are_similar__offset_1(const protein &prm_protein_a,         ///< TODOCUMENT
                                                 const protein &prm_protein_b,         ///< TODOCUMENT
                                                 const size_t  &prm_index_a__offset_1, ///< TODOCUMENT
                                                 const size_t  &prm_index_b__offset_1  ///< TODOCUMENT
                                                 ) const {
	check_offset_1( prm_index_a__offset_1 );
	check_offset_1( prm_index_b__offset_1 );

	const sec_struc &sec_struc_a = prm_protein_a.get_sec_struc_ref_of_index( prm_index_a__offset_1 - 1 );
	const sec_struc &sec_struc_b = prm_protein_b.get_sec_struc_ref_of_index( prm_index_b__offset_1 - 1 );

	return ( sec_struc_a.get_type() == sec_struc_b.get_type() );
}

/// \brief TODOCUMENT
bool sec_struc_querier::do_temp_hacky_is_residue() const {
	return false;
}


