/// \file
/// \brief Headers for convenience functions for handling pair alignments

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_PAIR_ALIGNMENT_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_PAIR_ALIGNMENT_HPP

#include "cath/alignment/alignment.hpp"

// clang-format off
namespace cath { class protein; }
namespace cath { class residue; }
// clang-format on

namespace cath::align {

	// Check that is a pair alignment
	void check_alignment_is_a_pair(const alignment &);

	// A bunch of convenience functions for pair alignments
	aln_posn_opt a_position_of_index(const alignment &, const size_t &);
	aln_posn_opt b_position_of_index(const alignment &, const size_t &);

	bool alignment_contains_pair(const alignment &, const size_t &, const size_t &);

	bool has_a_position_of_index(     const alignment &, const size_t & );
	bool has_b_position_of_index(     const alignment &, const size_t & );
	bool has_both_positions_of_index( const alignment &, const size_t & );

	aln_posn_type get_a_position_of_index(          const alignment &, const size_t & );
	aln_posn_type get_b_position_of_index(          const alignment &, const size_t & );
	aln_posn_type get_a_offset_1_position_of_index( const alignment &, const size_t & );
	aln_posn_type get_b_offset_1_position_of_index( const alignment &, const size_t & );

	bool has_last_a_position(        const alignment & );
	bool has_last_b_position(        const alignment & );
	bool has_last_position_of_entry( const alignment &, const size_t & );

	aln_posn_opt get_first_present_a_position(const alignment &);
	aln_posn_opt get_first_present_b_position(const alignment &);
	aln_posn_opt get_last_present_a_position(const alignment &);
	aln_posn_opt get_last_present_b_position(const alignment &);

	aln_posn_type get_last_a_position(                 const alignment & );
	aln_posn_type get_last_b_position(                 const alignment & );
	aln_posn_type get_last_position_of_entry(          const alignment &, const size_t & );
	aln_posn_type get_last_a_offset_1_position(        const alignment & );
	aln_posn_type get_last_b_offset_1_position(        const alignment & );
	aln_posn_type get_last_offset_1_position_of_entry( const alignment &, const size_t & );

	void append_position_a(             alignment &, const aln_posn_type & );
	void append_position_b(             alignment &, const aln_posn_type & );
	void append_position_both(          alignment &, const aln_posn_type &, const aln_posn_type & );
	void append_position_a_offset_1(    alignment &, const aln_posn_type & );
	void append_position_b_offset_1(    alignment &, const aln_posn_type & );
	void append_position_both_offset_1( alignment &, const aln_posn_type &, const aln_posn_type & );

	const residue & get_a_residue_cref_of_index(const alignment &, const protein &, const size_t &);
	const residue & get_b_residue_cref_of_index(const alignment &, const protein &, const size_t &);

	size_t num_present_positions_of_both_entries(const alignment &);
	size_vec indices_of_present_positions_of_both_entries(const alignment &);

	void set_pair_alignment_duplicate_scores(alignment &, const score_opt_vec &);
	alignment set_pair_alignment_duplicate_scores_copy(alignment, const score_opt_vec &);

} // namespace cath::align

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_PAIR_ALIGNMENT_HPP
