/// \file
/// \brief Definitions for convenience functions for handling pair alignments

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

#include "pair_alignment.h"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/copy_build.h"
#include "common/size_t_literal.h"
#include "common/temp_check_offset_1.h"
#include "exception/invalid_argument_exception.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::adaptors::filtered;
using boost::irange;

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::check_alignment_is_a_pair(const alignment &arg_alignment ///< TODOCUMENT
                                            ) {
	if (arg_alignment.num_entries() != alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Alignment is not a pair alignment"));
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment
opt_aln_posn cath::align::a_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                              const size_t    &arg_index      ///< The index to query
                                              ) {
	check_alignment_is_a_pair(arg_alignment);
	return arg_alignment.position_of_entry_of_index( alignment::PAIR_A_IDX, arg_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
opt_aln_posn cath::align::b_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                              const size_t    &arg_index      ///< The index to query
                                              ) {
	check_alignment_is_a_pair( arg_alignment );
	return arg_alignment.position_of_entry_of_index( alignment::PAIR_B_IDX, arg_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::alignment_contains_pair(const alignment &arg_alignment, ///< TODOCUMENT
                                          const size_t    &arg_posn_a,    ///< TODOCUMENT
							              const size_t    &arg_posn_b     ///< TODOCUMENT
							              ) {
	check_alignment_is_a_pair( arg_alignment );

	const size_t length = arg_alignment.length();
	for (size_t index = 0; index < length; ++index) {
		const opt_aln_posn &posn_a = arg_alignment.position_of_entry_of_index( alignment::PAIR_A_IDX, index );
		const opt_aln_posn &posn_b = arg_alignment.position_of_entry_of_index( alignment::PAIR_B_IDX, index );
		if ( posn_a && *posn_a == arg_posn_a && posn_b && *posn_b == arg_posn_b ) {
			return true;
		}
	}
	return false;
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_a_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                          const size_t    &arg_index      ///< The index to query
                                          ) {
	return static_cast<bool>( a_position_of_index( arg_alignment, arg_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_b_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                          const size_t    &arg_index      ///< The index to query
                                          ) {
	return static_cast<bool>( b_position_of_index( arg_alignment, arg_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_both_positions_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                              const size_t    &arg_index      ///< The index to query
                                              ) {
	return (
		has_a_position_of_index( arg_alignment, arg_index )
		&&
		has_b_position_of_index( arg_alignment, arg_index )
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_a_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                   const size_t    &arg_index      ///< The index to query
                                                   ) {
	// This uses get_position_of_entry_of_index() so that it benefits from the check that will
	// ensure errors result in an invalid_error_exception rather than an exception from boost::optional
	return get_position_of_entry_of_index( arg_alignment, alignment::PAIR_A_IDX, arg_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_b_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                   const size_t    &arg_index      ///< The index to query
                                                   ) {
	// This uses get_position_of_entry_of_index() so that it benefits from the check that will
	// ensure errors result in an invalid_error_exception rather than an exception from boost::optional
	return get_position_of_entry_of_index( arg_alignment, alignment::PAIR_B_IDX, arg_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
///
/// WARNING: The returned value starts at 1, not 0
aln_posn_type cath::align::get_a_offset_1_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                            const size_t    &arg_index      ///< The index to query
                                                            ) {
	return get_a_position_of_index(arg_alignment, arg_index) + 1;
}

/// \brief TODOCUMENT
///
/// \relates alignment
///
/// WARNING: The returned value starts at 1, not 0
aln_posn_type cath::align::get_b_offset_1_position_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                            const size_t    &arg_index      ///< The index to query
                                                            ) {
	return get_b_position_of_index(arg_alignment, arg_index) + 1;
}

///\brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_last_a_position(const alignment &arg_alignment ///< TODOCUMENT
                                      ) {
	return has_last_position_of_entry( arg_alignment, alignment::PAIR_A_IDX );
}

///\brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_last_b_position(const alignment &arg_alignment ///< TODOCUMENT
                                      ) {
	return has_last_position_of_entry( arg_alignment, alignment::PAIR_B_IDX );
}

///\brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_last_position_of_entry(const alignment &arg_alignment, ///< TODOCUMENT
                                             const size_t    &arg_entry      ///< TODOCUMENT
                                             ) {
	return has_position_of_entry_of_index( arg_alignment, arg_entry, arg_alignment.length() - 1);
}

/// \brief Convenience function for finding the index of the first present position for the first half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the first half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
opt_aln_posn cath::align::get_first_present_a_position(const alignment &arg_alignment ///< The alignment to be queried
                                                       ) {
	return get_first_present_position_of_entry( arg_alignment, alignment::PAIR_A_IDX );
}

/// \brief Convenience function for finding the index of the first present position for the second half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the second half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
opt_aln_posn cath::align::get_first_present_b_position(const alignment &arg_alignment ///< The alignment to be queried
                                                       ) {
	return get_first_present_position_of_entry( arg_alignment, alignment::PAIR_B_IDX );
}

/// \brief Convenience function for finding the index of the last present position for the first half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the first half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
opt_aln_posn cath::align::get_last_present_a_position(const alignment &arg_alignment ///< The alignment to be queried
                                                      ) {
	return get_last_present_position_of_entry( arg_alignment, alignment::PAIR_A_IDX );
}

/// \brief Convenience function for finding the index of the last present position for the second half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the second half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
opt_aln_posn cath::align::get_last_present_b_position(const alignment &arg_alignment ///< The alignment to be queried
                                                      ) {
	return get_last_present_position_of_entry( arg_alignment, alignment::PAIR_B_IDX );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_a_position(const alignment &arg_alignment ///< TODOCUMENT
                                               ) {
	return get_last_position_of_entry( arg_alignment, alignment::PAIR_A_IDX );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_b_position(const alignment &arg_alignment ///< TODOCUMENT
                                               ) {
	return get_last_position_of_entry( arg_alignment, alignment::PAIR_B_IDX );
}

///\brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_position_of_entry(const alignment &arg_alignment, ///< TODOCUMENT
                                                      const size_t    &arg_entry      ///< TODOCUMENT
                                                      ) {
	return get_position_of_entry_of_index( arg_alignment, arg_entry, arg_alignment.length() - 1 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_a_offset_1_position(const alignment &arg_alignment ///< TODOCUMENT
                                                        ) {
	return get_last_a_position( arg_alignment ) + 1;
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_b_offset_1_position(const alignment &arg_alignment ///< TODOCUMENT
                                                        ) {
	return get_last_b_position( arg_alignment ) + 1;
}

///\brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_offset_1_position_of_entry(const alignment &arg_alignment, ///< TODOCUMENT
                                                               const size_t    &arg_entry      ///< TODOCUMENT
                                                               ) {
	return get_last_position_of_entry( arg_alignment, arg_entry ) + 1;
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_a(alignment           &arg_alignment,   ///<TODOCUMENT
                                    const aln_posn_type &arg_position_val ///< TODOCUMENT
                                    ) {
	check_alignment_is_a_pair(arg_alignment);
	const size_t length = arg_alignment.length();
	arg_alignment.set_position_value(alignment::PAIR_A_IDX, length, arg_position_val);
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_b(alignment           &arg_alignment,   ///<TODOCUMENT
                                    const aln_posn_type &arg_position_val ///< TODOCUMENT
                                    ) {
	check_alignment_is_a_pair(arg_alignment);
	const size_t length = arg_alignment.length();
	arg_alignment.set_position_value(alignment::PAIR_B_IDX, length, arg_position_val);
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_both(alignment           &arg_alignment,      ///< TODOCUMENT
                                       const aln_posn_type &arg_a_position_val, ///< TODOCUMENT
                                       const aln_posn_type &arg_b_position_val  ///< TODOCUMENT
                                       ) {
	check_alignment_is_a_pair(arg_alignment);
	const size_t length = arg_alignment.length();
	arg_alignment.set_position_value(alignment::PAIR_A_IDX, length, arg_a_position_val);
	arg_alignment.set_position_value(alignment::PAIR_B_IDX, length, arg_b_position_val);
//	arg_alignment.append_position(
//		{       true       ,       true       },
//		{ position_a_value , position_b_value }
//	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_a_offset_1(alignment           &arg_alignment,   ///<TODOCUMENT
                                             const aln_posn_type &arg_position_val ///< TODOCUMENT
                                             ) {
	check_offset_1(arg_position_val);
	append_position_a( arg_alignment, arg_position_val - 1 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_b_offset_1(alignment           &arg_alignment,   ///<TODOCUMENT
                                             const aln_posn_type &arg_position_val ///< TODOCUMENT
                                             ) {
	check_offset_1(arg_position_val);
	append_position_b( arg_alignment, arg_position_val - 1 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_both_offset_1(alignment           &arg_alignment,      ///< TODOCUMENT
                                                const aln_posn_type &arg_a_position_val, ///< TODOCUMENT
                                                const aln_posn_type &arg_b_position_val  ///< TODOCUMENT
                                                ) {
	check_offset_1(arg_a_position_val);
	check_offset_1(arg_b_position_val);
	append_position_both( arg_alignment, arg_a_position_val - 1, arg_b_position_val - 1 );
}

/// \brief Convenience function for getting a const reference to an a residue using a pair alignment and index
///
/// \relates alignment
const residue & cath::align::get_a_residue_cref_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                         const protein   &arg_protein,   ///< TODOCUMENT
                                                         const size_t    &arg_index      ///< TODOCUMENT
                                                         ) {
	check_alignment_is_a_pair(arg_alignment);
	const aln_posn_type a_position__offset_1 = get_a_offset_1_position_of_index( arg_alignment, arg_index );
	check_offset_1(a_position__offset_1);
	return arg_protein.get_residue_ref_of_index(a_position__offset_1 - 1);
}

/// \brief Convenience function for getting a const reference to an a residue using a pair alignment and index
///
/// \relates alignment
const residue & cath::align::get_b_residue_cref_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                         const protein   &arg_protein,   ///< TODOCUMENT
                                                         const size_t    &arg_index      ///< TODOCUMENT
                                                         ) {
	check_alignment_is_a_pair(arg_alignment);
	const aln_posn_type b_position__offset_1 = get_b_offset_1_position_of_index( arg_alignment, arg_index );
	check_offset_1(b_position__offset_1);
	return arg_protein.get_residue_ref_of_index(b_position__offset_1 - 1);
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_t cath::align::num_present_positions_of_both_entries(const alignment &arg_alignment ///< TODOCUMENT
                                                          ) {
	check_alignment_is_a_pair(arg_alignment);
	return num_present_positions_of_both_entries(arg_alignment, alignment::PAIR_A_IDX, alignment::PAIR_B_IDX );
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_vec cath::align::indices_of_present_positions_of_both_entries(const alignment &arg_alignment ///< TODOCUMENT
                                                                   ) {
	return copy_build<size_vec>(
		irange( 0_z, arg_alignment.length() - 1 )
			| filtered( [&] (const size_t &x) { return has_both_positions_of_index( arg_alignment, x ); } )
	);
}


/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::set_pair_alignment_duplicate_scores(alignment           &arg_alignment, ///< TODOCUMENT
                                                      const opt_score_vec &arg_scores     ///< TODOCUMENT
                                                      ) {
	set_scores( arg_alignment, opt_score_vec_vec( alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT, arg_scores ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
alignment cath::align::set_pair_alignment_duplicate_scores_copy(alignment            arg_alignment, ///< TODOCUMENT
                                                                const opt_score_vec &arg_scores     ///< TODOCUMENT
                                                                ) {
	set_pair_alignment_duplicate_scores( arg_alignment, arg_scores );
	return arg_alignment;
}
