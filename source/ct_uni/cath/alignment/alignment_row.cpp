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

#include "alignment_row.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::std;

using ::boost::algorithm::any_of;
using ::boost::lexical_cast;

/// \brief Sanity check the specified entry value is within the range given the current number of entries
///        and throw an exception if not
///
/// \pre The specified entry value must be less than num_entries() else an invalid_argument_exception will be thrown
void alignment_row::sanity_check_entry(const size_t &prm_entry ///< The entry value to check
                                       ) const {
	if ( prm_entry >= num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Entry "
			+ lexical_cast<string>( prm_entry )
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
void alignment_row::check_has_position_of_entry(const size_t &prm_entry ///< The entry value to check
                                                ) const {
	if ( ! has_position_of_entry( *this, prm_entry) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot use position at entry "
			+ lexical_cast<string>( prm_entry )
			+ " of alignment_row because it does not have a position there"
		));
	}
}

/// \brief Ctor for alignment_row from a vector of optional alignment positions
alignment_row::alignment_row(aln_posn_opt_vec prm_positions ///< The vector of optional alignment positions from which to build the alignment row
                             ) : positions{ std::move( prm_positions ) } {
}

/// \brief Get the number of entries in this alignment row
size_t alignment_row::num_entries() const {
	return positions.size();
}

/// \brief Get the aln_posn_opt of the specified entry
aln_posn_opt alignment_row::position_of_entry(const size_t &prm_entry ///< The entry to query in the alignment_row
                                              ) const {
	sanity_check_entry( prm_entry );
	return positions[ prm_entry ];
}

/// \brief Convenience function to query whether a given position
///
/// \relates alignment_row
bool cath::align::has_position_of_entry(const alignment_row &prm_row,  ///< TODOCUMENT
                                        const size_t        &prm_entry ///< TODOCUMENT
                                        ) {
	return static_cast<bool>( prm_row.position_of_entry( prm_entry ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
aln_posn_type cath::align::get_position_of_entry(const alignment_row &prm_row,  ///< TODOCUMENT
                                                 const size_t        &prm_entry ///< TODOCUMENT
                                                 ) {
	return *prm_row.position_of_entry( prm_entry );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
///
/// \param prm_alignment_row TODOCUMENT
bool cath::align::any_entries_present( const alignment_row &prm_alignment_row ) {
	return any_of( indices( prm_alignment_row.num_entries() ),
	               [ & ]( const size_t &entry_ctr ) { return has_position_of_entry( prm_alignment_row, entry_ctr ); } );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
bool cath::align::any_entries_present(const alignment &prm_alignment, ///< TODOCUMENT
                                      const size_vec  &prm_entries,   ///< TODOCUMENT
                                      const size_t    &prm_index      ///< TODOCUMENT
                                      ) {
	return any_entries_present(
		get_row_of_entries_of_alignment(
			prm_alignment,
			prm_entries,
			prm_index
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::get_row_of_entries_of_alignment(const alignment &prm_alignment, ///< TODOCUMENT
                                                           const size_vec  &prm_entries,   ///< TODOCUMENT
                                                           const size_t    &prm_index      ///< TODOCUMENT
                                                           ) {
	const size_t num_entries = prm_entries.size();
	aln_posn_opt_vec positions;
	positions.reserve( num_entries );
	for (const size_t &entry : prm_entries) {
		const aln_posn_opt position = prm_alignment.position_of_entry_of_index( entry, prm_index );
		positions.push_back( position );
	}
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::get_row_of_alignment(const alignment &prm_alignment, ///< TODOCUMENT
                                                const size_t    &prm_index      ///< TODOCUMENT
                                                ) {
	const size_t num_entries = prm_alignment.num_entries();
	aln_posn_opt_vec positions;
	positions.reserve( num_entries );
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		const aln_posn_opt position = prm_alignment.position_of_entry_of_index( entry_ctr, prm_index );
		positions.push_back( position );
	}
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::make_empty_aln_row(const size_t &prm_num_empties ///< TODOCUMENT
                                              ) {
	return alignment_row( aln_posn_opt_vec ( prm_num_empties ) );
}

//		append_row( local_aln, make_row_with_single_value( entry, local_aln.num_entries, pdb_res_indices.size() - 1 ) );
//		add_row_with_single_value( local_aln, entry, pdb_res_indices.size() - 1 );

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::make_row_with_single_value(const size_t &prm_entry,       ///< TODOCUMENT
                                                      const size_t &prm_num_entries, ///< TODOCUMENT
                                                      const size_t &prm_value        ///< TODOCUMENT
                                                      ) {
	aln_posn_opt_vec positions( prm_num_entries );
	positions[ prm_entry ] = prm_value;
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
void cath::align::append_row_with_single_value(alignment    &prm_alignment, ///< TODOCUMENT
                                               const size_t &prm_entry,     ///< TODOCUMENT
                                               const size_t &prm_value      ///< TODOCUMENT
                                               ) {
	append_row(
		prm_alignment,
		make_row_with_single_value(
			prm_entry,
			prm_alignment.num_entries(),
			prm_value
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::copy_aln_row_entry(const alignment_row &prm_row,               ///< TODOCUMENT
                                              const size_t        &prm_copy_target_entry, ///< TODOCUMENT
                                              const alignment_row &prm_other_row,         ///< TODOCUMENT
                                              const size_t        &prm_copy_source_entry  ///< TODOCUMENT
                                              ) {
	aln_posn_opt_vec positions = get_has_posns_and_posns(prm_row);
	positions[ prm_copy_target_entry ] = prm_other_row.position_of_entry( prm_copy_source_entry );
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::glue_aln_row_onto_empties(const size_t        &prm_num_entries_in_empties,    ///< TODOCUMENT
                                                     const size_t        &prm_identify_entry_in_empties, ///< TODOCUMENT
                                                     const alignment_row &prm_aln_row,                   ///< TODOCUMENT
                                                     const size_t        &prm_identify_entry_in_row      ///< TODOCUMENT
                                                     ) {
	const alignment_row all_empties_row = make_empty_aln_row( prm_num_entries_in_empties );
	const alignment_row empties_row     = copy_aln_row_entry( all_empties_row, prm_identify_entry_in_empties, prm_aln_row, prm_identify_entry_in_row);
	const alignment_row row_from_b      = remove_entry_from_row( prm_aln_row, prm_identify_entry_in_row );
	return join( empties_row, row_from_b );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::glue_empties_onto_row(const alignment_row &prm_aln_row,               ///< TODOCUMENT
                                                 const size_t        &prm_num_entries_in_empties ///< TODOCUMENT
                                                 ) {
	return join( prm_aln_row, make_empty_aln_row( prm_num_entries_in_empties - 1 ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::glue_aln_rows_together(const alignment_row &prm_aln_row_a,        ///< TODOCUMENT
                                                  const size_t        &prm_identify_entry_a, ///< TODOCUMENT
                                                  const alignment_row &prm_aln_row_b,        ///< TODOCUMENT
                                                  const size_t        &prm_identify_entry_b  ///< TODOCUMENT
                                                  ) {
	// Grab the two positions and check they're both present
	const aln_posn_opt position_a = prm_aln_row_a.position_of_entry( prm_identify_entry_a );
	const aln_posn_opt position_b = prm_aln_row_b.position_of_entry( prm_identify_entry_b );
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
	return join( prm_aln_row_a, remove_entry_from_row( prm_aln_row_b, prm_identify_entry_b ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::join(const alignment_row &prm_joinee_a, ///< TODOCUMENT
                                const alignment_row &prm_joinee_b  ///< TODOCUMENT
                                ) {
	// Grab the total number of entries
	const size_t total_num_entries = prm_joinee_a.num_entries() + prm_joinee_b.num_entries();

	// Prepare a positions to be populated
	aln_posn_opt_vec positions;
	positions.reserve( total_num_entries );

	// Loop over the two joinees, populating positions
	for ( const alignment_row &joinee : { cref( prm_joinee_a ), cref( prm_joinee_b ) } ) {

		// Loop over the entries in the joinee
		const size_t num_entries = joinee.num_entries();
		for (const size_t &entry_ctr : indices( num_entries ) ) {
			// Append the position to positions
			const aln_posn_opt position = joinee.position_of_entry( entry_ctr );
			positions.push_back( position );
		}
	}

	// Return an alignment_row based on the positions
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
void cath::align::set_positions_of_entries_from_row(aln_posn_opt_vec    &prm_positions,  ///< TODOCUMENT
                                                    const alignment_row &prm_source_row, ///< TODOCUMENT
                                                    const size_vec      &prm_entries     ///< TODOCUMENT
                                                    ) {
	const size_t num_dest_entries   = prm_positions.size();
	const size_t num_source_entries = prm_entries.size();
	if ( prm_source_row.num_entries() != num_source_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of entries in does not match number of entries in entry list"));
	}

	for (const size_t &entry_ctr : indices( num_source_entries ) ) {
		const size_t new_entry = prm_entries[ entry_ctr ];
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
		if ( prm_positions[ new_entry ] ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to set a value that has already been set"));
		}
		prm_positions[ new_entry ] = prm_source_row.position_of_entry( entry_ctr );
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::weave(const alignment_row &prm_row_a,     ///< TODOCUMENT
                                 const size_vec      &prm_entries_a, ///< TODOCUMENT
                                 const alignment_row &prm_row_b,     ///< TODOCUMENT
                                 const size_vec      &prm_entries_b  ///< TODOCUMENT
                                 ) {
	const size_t num_entries_a = prm_entries_a.size();
	const size_t num_entries_b = prm_entries_b.size();
	const size_t num_entries = num_entries_a + num_entries_b;
	aln_posn_opt_vec positions( num_entries );
	set_positions_of_entries_from_row( positions, prm_row_a, prm_entries_a );
	set_positions_of_entries_from_row( positions, prm_row_b, prm_entries_b );
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
alignment_row cath::align::remove_entry_from_row( const alignment_row &prm_aln_row, ///< TODOCUMENT
                                                  const size_t &       prm_entry    ///< TODOCUMENT
                                                  ) {
	// Grab the number of entries and sanity check that prm_entry is less than this
	const size_t num_entries = prm_aln_row.num_entries();
	if ( prm_entry >= num_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot remove entry "
			+ lexical_cast<string>( prm_entry   )
			+ " from alignment_row with "
			+ lexical_cast<string>( num_entries )
			+ " entries"
		));
	}

	// Prepare a positions to be populated
	aln_posn_opt_vec positions;
	positions.reserve( num_entries );

	// Loop over the entries in prm_aln_row
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		// If this isn't the entry to be removed then copy it to has_positions and positions
		if ( prm_entry != entry_ctr ) {
			// Append the position to positions
			const aln_posn_opt position = prm_aln_row.position_of_entry( entry_ctr );
			positions.push_back( position );
		}
	}

	// Return an alignment_row based on the populated has_positions and positions
	return alignment_row( positions );
}

/// \brief TODOCUMENT
///
/// \relates alignment_row
aln_posn_opt_vec cath::align::get_has_posns_and_posns(const alignment_row &prm_aln_row ///< TODOCUMENT
                                                      ) {
	const size_t num_entries = prm_aln_row.num_entries();
	aln_posn_opt_vec new_positions( num_entries );
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		const aln_posn_opt position = prm_aln_row.position_of_entry( entry_ctr );
		new_positions[ entry_ctr ] = position;
	}
	return new_positions;
}

