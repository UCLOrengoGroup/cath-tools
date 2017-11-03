/// \file
/// \brief The alignment action definitions

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

#include "alignment_action.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>

#include "alignment/alignment.hpp"
#include "alignment/alignment_row.hpp"
#include "alignment/detail/multi_align_builder.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/exception/not_implemented_exception.hpp"

using namespace boost::log;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;

/// \brief A convenience function to append a glued row to the end of an alignment
void cath::align::detail::append_glued_row(alignment                  &arg_alignment,     ///< The alignment to which the glued row should be appended
                                           const aln_ent_ind_tup_pair &arg_data,          ///< The data structure containing the two source (alignment+entry+index)s
                                           const glued_row_type       &arg_glued_row_type ///< The type of glued row to add (from a, from b or from both)
                                           ) {
	// Grab the details from the aln_ent_ind_tup_pair
	const alignment &alignment_a = get<0>( arg_data.first  );
	const size_t    &entry_a     = get<1>( arg_data.first  );
	const size_t    &index_a     = get<2>( arg_data.first  );
	const alignment &alignment_b = get<0>( arg_data.second );
	const size_t    &entry_b     = get<1>( arg_data.second );
	const size_t    &index_b     = get<2>( arg_data.second );

	// Switch on the type of glue row (add a, add b or add both)
	switch (arg_glued_row_type) {
		case ( glued_row_type::FROM_A    ) : {
			append_row( arg_alignment, glue_empties_onto_row(
				get_row_of_alignment( alignment_a, index_a ),
				alignment_b.num_entries()
			));
			break;
		}
		case ( glued_row_type::FROM_B    ) : {
			append_row( arg_alignment, glue_aln_row_onto_empties(
				alignment_a.num_entries(),                           entry_a,
				get_row_of_alignment( alignment_b, index_b ), entry_b
			));
			break;
		}
		case ( glued_row_type::FROM_BOTH ) : {
			append_row( arg_alignment, glue_aln_rows_together(
				get_row_of_alignment( alignment_a, index_a ), entry_a,
				get_row_of_alignment( alignment_b, index_b ), entry_b
			));
			break;
		}
	}
}

/// \brief Glue two alignments by identifying one of the entries in the first alignment with one of the entries in the second alignment
///
/// \pre TODOCUMENT
///
/// \relates alignment
///
/// \returns TODOCUMENT
alignment cath::align::glue_two_alignments(const alignment &arg_alignment_a,    ///< The first alignment to be glued together
                                           const size_t    &arg_entry_in_aln_a, ///< The entry in the first  alignment to be identified with an entry in the second
                                           const alignment &arg_alignment_b,    ///< The second alignment to be glued together
                                           const size_t    &arg_entry_in_aln_b  ///< The entry in the second alignment to be identified with an entry in the first
                                           ) {
	using size_type = alignment::size_type;

	// Grab the number of entries and the lengths from the two alignments
	const size_type num_entries_in_a = arg_alignment_a.num_entries();
	const size_type num_entries_in_b = arg_alignment_b.num_entries();
	const size_type length_a         = arg_alignment_a.length();
	const size_type length_b         = arg_alignment_b.length();

	// Sanity check the indices of the two entries to be identified
	if ( arg_entry_in_aln_a >= num_entries_in_a || arg_entry_in_aln_b >= num_entries_in_b ) {
		const size_t &problem_index = ( arg_entry_in_aln_a >= num_entries_in_a ) ? arg_entry_in_aln_a   : arg_entry_in_aln_b;
		const size_t &problem_size  = ( arg_entry_in_aln_a >= num_entries_in_a ) ? num_entries_in_a : num_entries_in_b;
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to glue first alignment because index "
			+ lexical_cast<string>( problem_index )
			+ " is out of range in alignment with "
			+ lexical_cast<string>( problem_size  )
			+ " entries"
		));
	}

	// Create a new alignment with one fewer than the sum of the number of entries in the two alignments
	// (because two entries, one in each alignment, will be identified to make one new entry)
	alignment new_alignment( num_entries_in_a + num_entries_in_b - 1 );

	// Loop down the two alignments until the end of both has been reached,
	// keeping track of the position in the two alignments and the position in the entry being identified ("glue_ctr")
	size_type index_ctr_a = 0;
	size_type index_ctr_b = 0;
	while ( index_ctr_a < length_a || index_ctr_b < length_b ) {
		const bool before_end_of_a               = ( index_ctr_a < length_a );
		const bool before_end_of_b               = ( index_ctr_b < length_b );
		const bool alignment_a_has_glue_position = before_end_of_a && has_position_of_entry_of_index( arg_alignment_a, arg_entry_in_aln_a, index_ctr_a );
		const bool alignment_b_has_glue_position = before_end_of_b && has_position_of_entry_of_index( arg_alignment_b, arg_entry_in_aln_b, index_ctr_b );

		const aln_ent_ind_tup_pair a_and_b_data = make_pair(
			aln_ent_ind_tup( arg_alignment_a, arg_entry_in_aln_a, index_ctr_a ),
			aln_ent_ind_tup( arg_alignment_b, arg_entry_in_aln_b, index_ctr_b )
		);

		// If only the a side should be added, do that
		if ( before_end_of_a && ( ! before_end_of_b || ! alignment_a_has_glue_position ) ) {
			append_glued_row( new_alignment, a_and_b_data, glued_row_type::FROM_A );
			++index_ctr_a;
		}
		// If only the b side should be added, do that
		else if ( before_end_of_b && ( ! before_end_of_a || ! alignment_b_has_glue_position ) ) {
			append_glued_row( new_alignment, a_and_b_data, glued_row_type::FROM_B );
			++index_ctr_b;
		}
		// Otherwise both sides are before the end and both sides have a glue position so
		// both sides should be added - do that
		else {
			// Both positions should match so check that they do
			const aln_posn_type glue_position_a = get_position_of_entry_of_index( arg_alignment_a, arg_entry_in_aln_a, index_ctr_a );
			const aln_posn_type glue_position_b = get_position_of_entry_of_index( arg_alignment_b, arg_entry_in_aln_b, index_ctr_b );
			if ( glue_position_a == glue_position_b ) {

				// Both positions match so add both sides
				append_glued_row( new_alignment, a_and_b_data, glued_row_type::FROM_BOTH );
				++index_ctr_a;
				++index_ctr_b;
			}
			// Otherwise something funny is going on - a missing residue on one side of the alignment
			else {
				// Warn that problems have been detected
				BOOST_LOG_TRIVIAL( warning ) << "Whilst gluing alignments, found mismatching positions ( alignment_a[ entry: "
				                             << lexical_cast<string>( arg_entry_in_aln_a )
				                             << ", index: "
				                             << lexical_cast<string>( index_ctr_a        )
				                             << " ] at position "
				                             << lexical_cast<string>( glue_position_a    )
				                             << " and alignment_b[ entry: "
				                             << lexical_cast<string>( arg_entry_in_aln_b )
				                             << ", index: "
				                             << lexical_cast<string>( index_ctr_b        )
				                             << " ] at position "
				                             << lexical_cast<string>( glue_position_b    )
				                             << " )";


				// If glue_position_a is lower then add that side first
				if ( glue_position_a < glue_position_b ) {
					append_glued_row( new_alignment, a_and_b_data, glued_row_type::FROM_A );
					++index_ctr_a;
				}
				// If glue_position_b is lower then add that side first
				else {
					append_glued_row( new_alignment, a_and_b_data, glued_row_type::FROM_B );
					++index_ctr_b;
				}
			}
		}
	}

	return new_alignment;
}

/// \brief TODOCUMENT
///
/// \pre TODOCUMENT
///
/// \relates alignment
alignment cath::align::build_alignment_from_parts(const size_size_alignment_tuple_vec &arg_spanning_alignments, ///< TODOCUMENT
                                                  const protein_list                  &arg_proteins             ///< TODOCUMENT
                                                  ) {
	const size_t num_entries = arg_spanning_alignments.size() + 1;
	alignment new_alignment( num_entries );

	multi_align_builder builder( num_entries );

	for (const size_size_alignment_tuple &branch_alignment : arg_spanning_alignments) {
		add_alignment_branch( builder, branch_alignment, arg_proteins );
	}

	return builder.get_alignment();
}

