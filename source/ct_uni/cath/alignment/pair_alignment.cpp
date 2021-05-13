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

#include "pair_alignment.hpp"

#include <optional>

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/temp_check_offset_1.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;

using ::boost::adaptors::filtered;
using ::boost::algorithm::any_of;
using ::std::make_optional;

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::check_alignment_is_a_pair(const alignment &prm_alignment ///< TODOCUMENT
                                            ) {
	if (prm_alignment.num_entries() != alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Alignment is not a pair alignment"));
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_opt cath::align::a_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                              const size_t    &prm_index      ///< The index to query
                                              ) {
	check_alignment_is_a_pair(prm_alignment);
	return prm_alignment.position_of_entry_of_index( alignment::PAIR_A_IDX, prm_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_opt cath::align::b_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                              const size_t    &prm_index      ///< The index to query
                                              ) {
	check_alignment_is_a_pair( prm_alignment );
	return prm_alignment.position_of_entry_of_index( alignment::PAIR_B_IDX, prm_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
///
/// \param prm_alignment TODOCUMENT
/// \param prm_posn_a    TODOCUMENT
/// \param prm_posn_b    TODOCUMENT
bool cath::align::alignment_contains_pair( const alignment &prm_alignment, const size_t &prm_posn_a, const size_t &prm_posn_b ) {
	check_alignment_is_a_pair( prm_alignment );

	// clang-format off
	return any_of(
		indices( prm_alignment.length() ),
		[ & ]( const size_t &index ) {
			return (
				prm_alignment.position_of_entry_of_index( alignment::PAIR_A_IDX, index ) == make_optional( prm_posn_a )
				&&
				prm_alignment.position_of_entry_of_index( alignment::PAIR_B_IDX, index ) == make_optional( prm_posn_b )
			);
		}
	);
	// clang-format on
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_a_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                          const size_t    &prm_index      ///< The index to query
                                          ) {
	return static_cast<bool>( a_position_of_index( prm_alignment, prm_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_b_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                          const size_t    &prm_index      ///< The index to query
                                          ) {
	return static_cast<bool>( b_position_of_index( prm_alignment, prm_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_both_positions_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                              const size_t    &prm_index      ///< The index to query
                                              ) {
	return (
		has_a_position_of_index( prm_alignment, prm_index )
		&&
		has_b_position_of_index( prm_alignment, prm_index )
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_a_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                                   const size_t    &prm_index      ///< The index to query
                                                   ) {
	// This uses get_position_of_entry_of_index() so that it benefits from the check that will
	// ensure errors result in an invalid_error_exception rather than an exception from ::std::optional
	return get_position_of_entry_of_index( prm_alignment, alignment::PAIR_A_IDX, prm_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_b_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                                   const size_t    &prm_index      ///< The index to query
                                                   ) {
	// This uses get_position_of_entry_of_index() so that it benefits from the check that will
	// ensure errors result in an invalid_error_exception rather than an exception from ::std::optional
	return get_position_of_entry_of_index( prm_alignment, alignment::PAIR_B_IDX, prm_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
///
/// WARNING: The returned value starts at 1, not 0
aln_posn_type cath::align::get_a_offset_1_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                                            const size_t    &prm_index      ///< The index to query
                                                            ) {
	return get_a_position_of_index(prm_alignment, prm_index) + 1;
}

/// \brief TODOCUMENT
///
/// \relates alignment
///
/// WARNING: The returned value starts at 1, not 0
aln_posn_type cath::align::get_b_offset_1_position_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                                            const size_t    &prm_index      ///< The index to query
                                                            ) {
	return get_b_position_of_index(prm_alignment, prm_index) + 1;
}

///\brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_last_a_position(const alignment &prm_alignment ///< TODOCUMENT
                                      ) {
	return has_last_position_of_entry( prm_alignment, alignment::PAIR_A_IDX );
}

///\brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_last_b_position(const alignment &prm_alignment ///< TODOCUMENT
                                      ) {
	return has_last_position_of_entry( prm_alignment, alignment::PAIR_B_IDX );
}

///\brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_last_position_of_entry(const alignment &prm_alignment, ///< TODOCUMENT
                                             const size_t    &prm_entry      ///< TODOCUMENT
                                             ) {
	return has_position_of_entry_of_index( prm_alignment, prm_entry, prm_alignment.length() - 1);
}

/// \brief Convenience function for finding the index of the first present position for the first half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the first half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
aln_posn_opt cath::align::get_first_present_a_position(const alignment &prm_alignment ///< The alignment to be queried
                                                       ) {
	return get_first_present_position_of_entry( prm_alignment, alignment::PAIR_A_IDX );
}

/// \brief Convenience function for finding the index of the first present position for the second half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the second half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
aln_posn_opt cath::align::get_first_present_b_position(const alignment &prm_alignment ///< The alignment to be queried
                                                       ) {
	return get_first_present_position_of_entry( prm_alignment, alignment::PAIR_B_IDX );
}

/// \brief Convenience function for finding the index of the last present position for the first half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the first half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
aln_posn_opt cath::align::get_last_present_a_position(const alignment &prm_alignment ///< The alignment to be queried
                                                      ) {
	return get_last_present_position_of_entry( prm_alignment, alignment::PAIR_A_IDX );
}

/// \brief Convenience function for finding the index of the last present position for the second half in the specified alignment
///
/// \pre The alignment must contain at least one present position for the second half,
///      else an invalid_argument_exception will be thrown
///
/// \relates alignment
aln_posn_opt cath::align::get_last_present_b_position(const alignment &prm_alignment ///< The alignment to be queried
                                                      ) {
	return get_last_present_position_of_entry( prm_alignment, alignment::PAIR_B_IDX );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_a_position(const alignment &prm_alignment ///< TODOCUMENT
                                               ) {
	return get_last_position_of_entry( prm_alignment, alignment::PAIR_A_IDX );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_b_position(const alignment &prm_alignment ///< TODOCUMENT
                                               ) {
	return get_last_position_of_entry( prm_alignment, alignment::PAIR_B_IDX );
}

///\brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_position_of_entry(const alignment &prm_alignment, ///< TODOCUMENT
                                                      const size_t    &prm_entry      ///< TODOCUMENT
                                                      ) {
	return get_position_of_entry_of_index( prm_alignment, prm_entry, prm_alignment.length() - 1 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_a_offset_1_position(const alignment &prm_alignment ///< TODOCUMENT
                                                        ) {
	return get_last_a_position( prm_alignment ) + 1;
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_b_offset_1_position(const alignment &prm_alignment ///< TODOCUMENT
                                                        ) {
	return get_last_b_position( prm_alignment ) + 1;
}

///\brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_last_offset_1_position_of_entry(const alignment &prm_alignment, ///< TODOCUMENT
                                                               const size_t    &prm_entry      ///< TODOCUMENT
                                                               ) {
	return get_last_position_of_entry( prm_alignment, prm_entry ) + 1;
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_a(alignment           &prm_alignment,   ///<TODOCUMENT
                                    const aln_posn_type &prm_position_val ///< TODOCUMENT
                                    ) {
	check_alignment_is_a_pair(prm_alignment);
	const size_t length = prm_alignment.length();
	prm_alignment.set_position_value(alignment::PAIR_A_IDX, length, prm_position_val);
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_b(alignment           &prm_alignment,   ///<TODOCUMENT
                                    const aln_posn_type &prm_position_val ///< TODOCUMENT
                                    ) {
	check_alignment_is_a_pair(prm_alignment);
	const size_t length = prm_alignment.length();
	prm_alignment.set_position_value(alignment::PAIR_B_IDX, length, prm_position_val);
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_both(alignment           &prm_alignment,      ///< TODOCUMENT
                                       const aln_posn_type &prm_a_position_val, ///< TODOCUMENT
                                       const aln_posn_type &prm_b_position_val  ///< TODOCUMENT
                                       ) {
	check_alignment_is_a_pair(prm_alignment);
	const size_t length = prm_alignment.length();
	prm_alignment.set_position_value(alignment::PAIR_A_IDX, length, prm_a_position_val);
	prm_alignment.set_position_value(alignment::PAIR_B_IDX, length, prm_b_position_val);
//	prm_alignment.append_position(
//		{       true       ,       true       },
//		{ position_a_value , position_b_value }
//	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_a_offset_1(alignment           &prm_alignment,   ///<TODOCUMENT
                                             const aln_posn_type &prm_position_val ///< TODOCUMENT
                                             ) {
	check_offset_1(prm_position_val);
	append_position_a( prm_alignment, prm_position_val - 1 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_b_offset_1(alignment           &prm_alignment,   ///<TODOCUMENT
                                             const aln_posn_type &prm_position_val ///< TODOCUMENT
                                             ) {
	check_offset_1(prm_position_val);
	append_position_b( prm_alignment, prm_position_val - 1 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_position_both_offset_1(alignment           &prm_alignment,      ///< TODOCUMENT
                                                const aln_posn_type &prm_a_position_val, ///< TODOCUMENT
                                                const aln_posn_type &prm_b_position_val  ///< TODOCUMENT
                                                ) {
	check_offset_1(prm_a_position_val);
	check_offset_1(prm_b_position_val);
	append_position_both( prm_alignment, prm_a_position_val - 1, prm_b_position_val - 1 );
}

/// \brief Convenience function for getting a const reference to an a residue using a pair alignment and index
///
/// \relates alignment
const residue & cath::align::get_a_residue_cref_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                                         const protein   &prm_protein,   ///< TODOCUMENT
                                                         const size_t    &prm_index      ///< TODOCUMENT
                                                         ) {
	check_alignment_is_a_pair(prm_alignment);
	const aln_posn_type a_position__offset_1 = get_a_offset_1_position_of_index( prm_alignment, prm_index );
	check_offset_1(a_position__offset_1);
	return prm_protein.get_residue_ref_of_index(a_position__offset_1 - 1);
}

/// \brief Convenience function for getting a const reference to an a residue using a pair alignment and index
///
/// \relates alignment
const residue & cath::align::get_b_residue_cref_of_index(const alignment &prm_alignment, ///< TODOCUMENT
                                                         const protein   &prm_protein,   ///< TODOCUMENT
                                                         const size_t    &prm_index      ///< TODOCUMENT
                                                         ) {
	check_alignment_is_a_pair(prm_alignment);
	const aln_posn_type b_position__offset_1 = get_b_offset_1_position_of_index( prm_alignment, prm_index );
	check_offset_1(b_position__offset_1);
	return prm_protein.get_residue_ref_of_index(b_position__offset_1 - 1);
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_t cath::align::num_present_positions_of_both_entries(const alignment &prm_alignment ///< TODOCUMENT
                                                          ) {
	check_alignment_is_a_pair(prm_alignment);
	return num_present_positions_of_both_entries(prm_alignment, alignment::PAIR_A_IDX, alignment::PAIR_B_IDX );
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_vec cath::align::indices_of_present_positions_of_both_entries(const alignment &prm_alignment ///< TODOCUMENT
                                                                   ) {
	return copy_build<size_vec>(
		indices( prm_alignment.length() - 1 )
			| filtered( [&] (const size_t &x) { return has_both_positions_of_index( prm_alignment, x ); } )
	);
}


/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::set_pair_alignment_duplicate_scores(alignment           &prm_alignment, ///< TODOCUMENT
                                                      const score_opt_vec &prm_scores     ///< TODOCUMENT
                                                      ) {
	set_scores( prm_alignment, score_opt_vec_vec( alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT, prm_scores ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
alignment cath::align::set_pair_alignment_duplicate_scores_copy(alignment            prm_alignment, ///< TODOCUMENT
                                                                const score_opt_vec &prm_scores     ///< TODOCUMENT
                                                                ) {
	set_pair_alignment_duplicate_scores( prm_alignment, prm_scores );
	return prm_alignment;
}
