/// \file
/// \brief The alignment_row class definitions

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

#include "alignment_row.h"

#include <boost/lexical_cast.hpp>

#include "alignment/alignment.h"
#include "exception/invalid_argument_exception.h"

using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;

/// \brief Sanity check the specified entry value is within the range given the current number of entries
///        and throw an exception if not
///
/// \pre The specified entry value must be less than num_entries() else an invalid_argument_exception will be thrown
void alignment_row::sanity_check_entry(const size_t &arg_entry ///< The entry value to check
                                       ) const {
	if ( arg_entry >= num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Entry "
			+ lexical_cast<string>( arg_entry )
			+ " is out of range in this alignment_row, which has "
			+ lexical_cast<string>( num_entries() )
			+ " entries"
		));
	}
}

/// \brief Check that this alignment row has a position present at the specified entry value
///        and throw an exception if not
///
/// \pre This alignment row must has_position_of_entry() at this entry else an invalid_argument_exception will be thrown
void alignment_row::check_has_position_of_entry(const size_t &arg_entry ///< The entry value to check
                                                ) const {
	if ( ! has_position_of_entry( *this, arg_entry) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot use position at entry "
			+ lexical_cast<string>( arg_entry )
			+ " of alignment_row because it does not have a position there"
		));
	}
}

/// \brief Ctor for alignment_row from a vector of optional alignment positions
alignment_row::alignment_row(const opt_aln_posn_vec &arg_positions ///< The vector of optional alignment positions from which to build the alignment row
                             ) : positions ( arg_positions     ) {
}

/// \brief Get the number of entries in this alignment row
size_t alignment_row::num_entries() const {
	return positions.size();
}

/// \brief Get the opt_aln_posn of the specified entry
opt_aln_posn alignment_row::position_of_entry(const size_t &arg_entry ///< The entry to query in the alignment_row
                                              ) const {
	sanity_check_entry( arg_entry );
	return positions[ arg_entry ];
}

/// \brief Convenience function to query whether a given position
///
/// \relates alignment_row
bool cath::align::has_position_of_entry(const alignment_row &arg_row,  ///< TODOCUMENT
                                        const size_t        &arg_entry ///< TODOCUMENT
                                        ) {
	return static_cast<bool>( arg_row.position_of_entry( arg_entry ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
aln_posn_type cath::align::get_position_of_entry(const alignment_row &arg_row,  ///< TODOCUMENT
                                                 const size_t        &arg_entry ///< TODOCUMENT
                                                 ) {
	return *arg_row.position_of_entry( arg_entry );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
bool cath::align::any_entries_present(const alignment_row &arg_alignment_row ///< TODOCUMENT
                                      ) {
	const size_t num_entries = arg_alignment_row.num_entries();
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		if ( has_position_of_entry( arg_alignment_row, entry_ctr ) ) {
			return true;
		}
	}
	return false;
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
bool cath::align::any_entries_present(const alignment &arg_alignment, ///< TODOCUMENT
                                      const size_vec  &arg_entries,   ///< TODOCUMENT
                                      const size_t    &arg_index      ///< TODOCUMENT
                                      ) {
	return any_entries_present(
		get_row_of_entries_of_alignment(
			arg_alignment,
			arg_entries,
			arg_index
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::get_row_of_entries_of_alignment(const alignment &arg_alignment, ///< TODOCUMENT
                                                           const size_vec  &arg_entries,   ///< TODOCUMENT
                                                           const size_t    &arg_index      ///< TODOCUMENT
                                                           ) {
	const size_t num_entries = arg_entries.size();
	opt_aln_posn_vec positions;
	positions.reserve( num_entries );
	for (const size_t &entry : arg_entries) {
		const opt_aln_posn position = arg_alignment.position_of_entry_of_index( entry, arg_index );
		positions.push_back( position );
	}
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::get_row_of_alignment(const alignment &arg_alignment, ///< TODOCUMENT
                                                const size_t    &arg_index      ///< TODOCUMENT
                                                ) {
	const size_t num_entries = arg_alignment.num_entries();
	opt_aln_posn_vec positions;
	positions.reserve( num_entries );
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		const opt_aln_posn position = arg_alignment.position_of_entry_of_index( entry_ctr, arg_index );
		positions.push_back( position );
	}
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::make_empty_aln_row(const size_t &arg_num_empties ///< TODOCUMENT
                                              ) {
	return alignment_row( opt_aln_posn_vec ( arg_num_empties ) );
}

//		append_row( local_aln, make_row_with_single_value( entry, local_aln.num_entries, pdb_res_indices.size() - 1 ) );
//		add_row_with_single_value( local_aln, entry, pdb_res_indices.size() - 1 );

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::make_row_with_single_value(const size_t &arg_entry,       ///< TODOCUMENT
                                                      const size_t &arg_num_entries, ///< TODOCUMENT
                                                      const size_t &arg_value        ///< TODOCUMENT
                                                      ) {
	opt_aln_posn_vec positions( arg_num_entries );
	positions[ arg_entry ] = arg_value;
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
void cath::align::append_row_with_single_value(alignment    &arg_alignment, ///< TODOCUMENT
                                               const size_t &arg_entry,     ///< TODOCUMENT
                                               const size_t &arg_value      ///< TODOCUMENT
                                               ) {
	append_row(
		arg_alignment,
		make_row_with_single_value(
			arg_entry,
			arg_alignment.num_entries(),
			arg_value
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::copy_aln_row_entry(const alignment_row &arg_row,               ///< TODOCUMENT
                                              const size_t        &arg_copy_target_entry, ///< TODOCUMENT
                                              const alignment_row &arg_other_row,         ///< TODOCUMENT
                                              const size_t        &arg_copy_source_entry  ///< TODOCUMENT
                                              ) {
	opt_aln_posn_vec positions = get_has_posns_and_posns(arg_row);
	positions[ arg_copy_target_entry ] = arg_other_row.position_of_entry( arg_copy_source_entry );
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::glue_aln_row_onto_empties(const size_t        &arg_num_entries_in_empties,    ///< TODOCUMENT
                                                     const size_t        &arg_identify_entry_in_empties, ///< TODOCUMENT
                                                     const alignment_row &arg_aln_row,                   ///< TODOCUMENT
                                                     const size_t        &arg_identify_entry_in_row      ///< TODOCUMENT
                                                     ) {
	const alignment_row all_empties_row = make_empty_aln_row( arg_num_entries_in_empties );
	const alignment_row empties_row     = copy_aln_row_entry( all_empties_row, arg_identify_entry_in_empties, arg_aln_row, arg_identify_entry_in_row);
	const alignment_row row_from_b      = remove_entry_from_row( arg_aln_row, arg_identify_entry_in_row );
	return join( empties_row, row_from_b );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::glue_empties_onto_row(const alignment_row &arg_aln_row,               ///< TODOCUMENT
                                                 const size_t        &arg_num_entries_in_empties ///< TODOCUMENT
                                                 ) {
	return join( arg_aln_row, make_empty_aln_row( arg_num_entries_in_empties - 1 ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::glue_aln_rows_together(const alignment_row &arg_aln_row_a,        ///< TODOCUMENT
                                                  const size_t        &arg_identify_entry_a, ///< TODOCUMENT
                                                  const alignment_row &arg_aln_row_b,        ///< TODOCUMENT
                                                  const size_t        &arg_identify_entry_b  ///< TODOCUMENT
                                                  ) {
	// Grab the two positions and check they're both present
	const opt_aln_posn position_a = arg_aln_row_a.position_of_entry( arg_identify_entry_a );
	const opt_aln_posn position_b = arg_aln_row_b.position_of_entry( arg_identify_entry_b );
	if ( ! position_a ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("first alignment_row does not have a position in the entry on which it is to be glued"));
	}
	if ( ! position_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("second alignment_row does not have a position in the entry on which it is to be glued"));
	}

	// Check the two positions both match
	if ( *position_a != *position_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"The positions on which two alignment_rows are to be glued do not match ( "
			+ lexical_cast<string>( *position_a )
			+ " != "
			+ lexical_cast<string>( *position_b )
			+ " )"
		));
	}

	// Return the result
	return join( arg_aln_row_a, remove_entry_from_row( arg_aln_row_b, arg_identify_entry_b ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::join(const alignment_row &arg_joinee_a, ///< TODOCUMENT
                                const alignment_row &arg_joinee_b  ///< TODOCUMENT
                                ) {
	// Grab the total number of entries
	const size_t num_entries = arg_joinee_a.num_entries() + arg_joinee_b.num_entries();

	// Prepare a positions to be populated
	opt_aln_posn_vec positions;
	positions.reserve( num_entries );

	// Loop over the two joinees, populating positions
	for ( const alignment_row &joinee : { cref( arg_joinee_a ), cref( arg_joinee_b ) } ) {

		// Loop over the entries in the joinee
		const size_t num_entries = joinee.num_entries();
		for ( size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr ) {
			// Append the position to positions
			const opt_aln_posn position = joinee.position_of_entry( entry_ctr );
			positions.push_back( position );
		}
	}

	// Return a alignment_row based on the positions
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
void cath::align::set_positions_of_entries_from_row(opt_aln_posn_vec    &arg_positions,  ///< TODOCUMENT
                                                    const alignment_row &arg_source_row, ///< TODOCUMENT
                                                    const size_vec      &arg_entries     ///< TODOCUMENT
                                                    ) {
	const size_t num_dest_entries   = arg_positions.size();
	const size_t num_source_entries = arg_entries.size();
	if ( arg_source_row.num_entries() != num_source_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of entries in does not match number of entries in entry list"));
	}

	for (size_t entry_ctr = 0; entry_ctr < num_source_entries; ++entry_ctr) {
		const size_t new_entry = arg_entries[ entry_ctr ];
//		cerr << "entry_ctr is " << entry_ctr << " and new_entry is " << new_entry << endl;
		if ( new_entry >= num_dest_entries ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"When trying to set position with number "
				+ lexical_cast<string>( entry_ctr )
				+ " in the source, the entry number in the destination is "
				+ lexical_cast<string>( new_entry )
				+ ", which is out of range because there are only "
				+ lexical_cast<string>( num_dest_entries )
				+ " destination entries"
			));
		}
		if ( arg_positions[ new_entry ] ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to set a value that has already been set"));
		}
		arg_positions[ new_entry ] = arg_source_row.position_of_entry( entry_ctr );
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::weave(const alignment_row &arg_row_a,     ///< TODOCUMENT
                                 const size_vec      &arg_entries_a, ///< TODOCUMENT
                                 const alignment_row &arg_row_b,     ///< TODOCUMENT
                                 const size_vec      &arg_entries_b  ///< TODOCUMENT
                                 ) {
	const size_t num_entries_a = arg_entries_a.size();
	const size_t num_entries_b = arg_entries_b.size();
	const size_t num_entries = num_entries_a + num_entries_b;
	opt_aln_posn_vec positions( num_entries );
	set_positions_of_entries_from_row( positions, arg_row_a, arg_entries_a );
	set_positions_of_entries_from_row( positions, arg_row_b, arg_entries_b );
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
const alignment_row cath::align::remove_entry_from_row(const alignment_row &arg_aln_row, ///< TODOCUMENT
                                                       const size_t        &arg_entry    ///< TODOCUMENT
                                                       ) {
	// Grab the number of entries and sanity check that arg_entry is less than this
	const size_t num_entries = arg_aln_row.num_entries();
	if ( arg_entry >= num_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot remove entry "
			+ lexical_cast<string>( arg_entry   )
			+ " from alignment_row with "
			+ lexical_cast<string>( num_entries )
			+ " entries"
		));
	}

	// Prepare a positions to be populated
	opt_aln_posn_vec positions;
	positions.reserve( num_entries );

	// Loop over the entries in arg_aln_row
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		// If this isn't the entry to be removed then copy it to has_positions and positions
		if ( arg_entry != entry_ctr ) {
			// Append the position to positions
			const opt_aln_posn position = arg_aln_row.position_of_entry( entry_ctr );
			positions.push_back( position );
		}
	}

	// Return a alignment_row based on the populated has_positions and positions
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
opt_aln_posn_vec cath::align::get_has_posns_and_posns(const alignment_row &arg_aln_row ///< TODOCUMENT
                                                      ) {
	const size_t num_entries = arg_aln_row.num_entries();
	opt_aln_posn_vec new_positions( num_entries );
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		const opt_aln_posn position = arg_aln_row.position_of_entry( entry_ctr );
		new_positions[ entry_ctr ] = position;
	}
	return new_positions;
}

