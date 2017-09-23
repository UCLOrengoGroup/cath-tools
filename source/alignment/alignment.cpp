/// \file
/// \brief The alignment class definitions

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

#include "alignment.hpp"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>
#include <boost/throw_exception.hpp>

#include "alignment/alignment_row.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/size_t_literal.hpp"
#include "common/temp_check_offset_1.hpp"
#include "exception/invalid_argument_exception.hpp"

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::all_of;
using boost::algorithm::any_of;
using boost::accumulate;
using boost::irange;
using boost::lexical_cast;
using boost::none;
using boost::numeric_cast;
using boost::range::count_if;
using boost::range::max_element;

constexpr alignment::size_type alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT;
constexpr alignment::size_type alignment::PAIR_A_IDX;
constexpr alignment::size_type alignment::PAIR_B_IDX;

/// \brief TODOCUMENT
void alignment::check_scored() const {
	if (!is_scored()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Alignment has not been scored"));
	}
}

/// \brief TODOCUMENT
void alignment::check_not_scored() const {
	if (is_scored()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to modify an alignment that has been scored"));
	}
}

/// \brief TODOCUMENT
void alignment::check_entry_in_range(const alignment::size_type &arg_entry ///< TODOCUMENT
                                     ) const {
	if (arg_entry >= num_entries()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Entry argument is greater than or equal to the number of entries being aligned"));
	}
}

/// \brief TODOCUMENT
void alignment::check_index_in_range(const alignment::size_type &arg_index ///< TODOCUMENT
                                     ) const {
	if (arg_index >= length()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Index argument is greater than or equal to the length of the alignment"));
	}
}

/// \brief TODOCUMENT
alignment::size_type alignment::reserved_length() const {
	return positions.empty() ? 0 : positions.front().size();
}

/// \brief Ctor for alignment
///
/// It is useful to allow an alignment with one entry so that cath-superpose doesn't
/// have to treat that as a special case.
alignment::alignment(const size_type &arg_num_entries ///< TODOCUMENT
                     ) : positions      ( arg_num_entries ),
                         logical_length ( 0               ) {
	if ( arg_num_entries < 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot currently create alignment of with no entries"));
	}
}

/// \brief A constructor to make an alignment from vectors of bool and aln_posn_type
///
/// It is useful to allow an alignment with one entry so that cath-superpose doesn't
/// have to treat that as a special case.
alignment::alignment(const aln_posn_opt_vec_vec &arg_lists ///< TODOCUMENT
                     ) : positions      ( arg_lists.size() ),
                         logical_length ( 0                ) {
	// If there isn't at least one list, then throw a wobbly
	const size_type the_num_entries = arg_lists.size();
	if (the_num_entries < 1) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct an alignment from zero lists"));
	}

	// If any of the lists are of different sizes, then throw a wobbly
	const size_type input_length = arg_lists.front().size();
	for (const aln_posn_opt_vec &list : arg_lists) {
		if ( input_length != list.size() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct alignment from lists of differing lengths"));
		}
	}

	// Reserve the space for these positions
	reserve( input_length );

	// Step through the alignment, adding the new bunch of positions
	for (size_t index = 0; index < input_length; ++index ) {
		aln_posn_opt_vec new_positions( the_num_entries );

		bool has_position_for_some_entry = false;
		for (alignment::size_type entry = 0; entry < the_num_entries; ++entry) {
			const aln_posn_opt &position = arg_lists[ entry ][ index ];
			has_position_for_some_entry  = ( position || has_position_for_some_entry );
			new_positions[ entry ]       = position;
		}

		if ( ! has_position_for_some_entry )  {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct an alignment since there is an index for which no entry has a position"));
		}

//		cerr << "About to append position for " << aln_ctr << endl;
		append_row( *this, new_positions );
	}
}

/// \brief TODOCUMENT
void alignment::reserve(const size_type &arg_size ///< TODOCUMENT
                        ) {
	for (aln_posn_opt_vec &position_part : positions) {
		position_part.resize( arg_size, none );
	}
}

/// \brief TODOCUMENT
alignment::size_type alignment::num_entries() const {
	return positions.size();
}

/// \brief TODOCUMENT
alignment::size_type alignment::length() const {
	return logical_length;
}

/// \brief TODOCUMENT
bool alignment::is_scored() const {
	return static_cast<bool>( new_scores );
}

/// \brief TODOCUMENT
const alignment_residue_scores & alignment::get_alignment_residue_scores() const {
	return *new_scores;
}

/// \brief TODOCUMENT
aln_posn_opt alignment::position_of_entry_of_index(const size_type &arg_entry, ///< TODOCUMENT
                                                   const size_type &arg_index  ///< TODOCUMENT
                                                   ) const {
	check_entry_in_range( arg_entry );
	check_index_in_range( arg_index );
	return positions[ arg_entry ][ arg_index ];
}

/// \brief TODOCUMENT
void alignment::set_position_value(const size_type     &arg_entry, ///< TODOCUMENT
                                   const size_type     &arg_index, ///< TODOCUMENT
                                   const aln_posn_type &arg_value  ///< TODOCUMENT
                                   ) {
	// Check that this alignment hasn't already been scored (in which case it should be read-only)
	check_not_scored();

	const size_type aln_num_entries = num_entries();

	if (arg_entry >= aln_num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Not currently able to add new entries to existing alignment"));
	}
	if (arg_index > length()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Not currently able to add multiple rows to an existing alignment"));
	}

	if ( arg_index >= reserved_length() ) {
		for (aln_posn_opt_vec &position_part : positions ) {
			position_part.resize ( arg_index + 1, none );
		}
	}

	if (arg_index >= length()) {
		logical_length = arg_index+1;
	}

	// Set the entry
	positions[ arg_entry ][ arg_index ] = arg_value;
}

/// \brief TODOCUMENT
void alignment::set_scores(const alignment_residue_scores &arg_scores ///< TODOCUMENT
                           ) {
	if ( arg_scores.get_num_entries() != num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of entries in scores does not match the number of entries in the alignment"));
	}
	if ( arg_scores.get_length() != length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Length of scores does not match the length of the alignment"));
	}
	const size_t the_length = length();
	for (size_t index_ctr = 0; index_ctr < the_length; ++index_ctr) {
		if ( arg_scores.get_num_present_entries_of_index( index_ctr ) != num_present_positions_of_index( *this, index_ctr ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of present positions at index does not match the equivalent in the alignment"));
		}
	}
	new_scores = arg_scores;
}

/// \brief Return the first non-consecutive position in any entry of the alignment, if any
///
/// To be clearer about "first": all positions in one entry are considered to be before any in the next entry
///
/// \returns A pair of the entry and expected position of the first non-consecutive entry, if any,
///          or none otherwise.
///
/// \relates alignment
size_size_pair_opt cath::align::first_non_consecutive_entry_positions(const alignment &arg_alignment ///< The alignment to search
                                                                      ) {
	const alignment::size_type num_entries = arg_alignment.num_entries();
	const alignment::size_type length      = arg_alignment.length();

	// Loop over the entries
	for (size_t entry = 0; entry < num_entries; ++entry) {

		// Loop over the indices for the entry, keeping track of the expected positions
		size_t expected_position = 0;
		for (size_t index = 0; index < length; ++index) {
			// If there is a position...
			if ( has_position_of_entry_of_index( arg_alignment, entry,index ) ) {
				// If the found position doesn't match the expected position, return the entry and expected position
				const aln_posn_type got_position = get_position_of_entry_of_index( arg_alignment, entry,index );
				if ( got_position != expected_position ) {
					return make_pair( entry, expected_position );
				}

				// Increment the expected position
				++expected_position;
			}
		}
	}

	// No non-consecutive position was found, so return none
	return none;
}

/// \brief Check that alignments positions are consecutive for each entry, else throw an exception
///
/// \pre The alignment may not contain any non-consecutive positions, else an invalid_argument_exception_will be thrown
///
/// \relates alignment
void cath::align::check_entry_positions_are_consecutive(const alignment &arg_alignment ///< The alignment to check
                                                        ) {
	// Grab the first non-consecutive position, if any
	const size_size_pair_opt non_consecutive = first_non_consecutive_entry_positions( arg_alignment );

	// If any were found, then throw a descriptive exception
	if ( non_consecutive ) {
		const size_t &non_consecutive_entry    = non_consecutive->first;
		const size_t &non_consecutive_position = non_consecutive->second;
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Alignment has non-consecutive elements (in entry "
			+ lexical_cast<string>( non_consecutive_entry )
			+ " where position "
			+ lexical_cast<string>( non_consecutive_position )
			+ " was expected)"
		));
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_position_of_entry_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                 const size_t    &arg_entry,     ///< TODOCUMENT
                                                 const size_t    &arg_index      ///< TODOCUMENT
                                                 ) {
	return static_cast<bool>( arg_alignment.position_of_entry_of_index( arg_entry, arg_index ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_position_of_both_entries_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                        const size_t    &arg_entry_a,   ///< TODOCUMENT
                                                        const size_t    &arg_entry_b,   ///< TODOCUMENT
                                                        const size_t    &arg_index      ///< TODOCUMENT
                                                        ) {
	return (
		has_position_of_entry_of_index( arg_alignment, arg_entry_a, arg_index )
		&&
		has_position_of_entry_of_index( arg_alignment, arg_entry_b, arg_index )
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
bool cath::align::has_position_of_all_entries_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                       const size_t    &arg_index      ///< TODOCUMENT
                                                       ) {
	return all_of(
		irange( 0_z, arg_alignment.num_entries() ),
		[&] (const size_t &x) { return has_position_of_entry_of_index( arg_alignment, x, arg_index ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_type cath::align::get_position_of_entry_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                          const size_t    &arg_entry,     ///< TODOCUMENT
                                                          const size_t    &arg_index      ///< TODOCUMENT
                                                          ) {
	const aln_posn_opt position = arg_alignment.position_of_entry_of_index( arg_entry, arg_index );
	if ( ! position ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("No position value present in entry " + lexical_cast<string>(arg_entry) + " at index " + lexical_cast<string>(arg_index)));
	}
	return *position;
}

/// \brief Whether there are any present positions in the specified range of consecutive indices in
///        the specified entry in the specified alignment
///
/// \relates alignment
bool cath::align::has_positions_of_entry_in_index_range(const alignment &arg_alignment,   ///< The alignment to be inspected
                                                        const size_t    &arg_entry,       ///< The entry to be inspected
                                                        const size_t    &arg_begin_index, ///< Begin index of the range to be inspected
                                                        const size_t    &arg_end_index    ///< One-past-end index of the range to be inspected
                                                        ) {
	if ( arg_begin_index > arg_end_index ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot query an invalid range with begin > end"));
	}
	return any_of(
		irange( arg_begin_index, arg_end_index ),
		[&] (const size_t &x) {
			return has_position_of_entry_of_index( arg_alignment, arg_entry, x );
		}
	);
}

/// \brief Non-member equality operator for the alignment class
///
/// This does not current compare (the existence of scores) so this may return true for
/// two alignments where only one has scores or their scores differ
///
/// \relates alignment
bool cath::align::operator==(const alignment &arg_aln_a, ///< TODOCUMENT
                             const alignment &arg_aln_b  ///< TODOCUMENT
                             ) {
	// Return false if the lengths or number of entries differ
	if (arg_aln_a.length() != arg_aln_b.length()) {
		return false;
	}
	if (arg_aln_a.num_entries() != arg_aln_b.num_entries()) {
		return false;
	}

	// Otherwise, check each of the positions match
	const alignment::size_type length      = arg_aln_a.length();
	const alignment::size_type num_entries = arg_aln_a.num_entries();
	for (alignment::size_type index_ctr = 0; index_ctr < length; ++index_ctr) {
		for (alignment::size_type entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
			const aln_posn_opt &a_position = arg_aln_a.position_of_entry_of_index( entry_ctr, index_ctr );
			const aln_posn_opt &b_position = arg_aln_b.position_of_entry_of_index( entry_ctr, index_ctr );
			if ( a_position != b_position ) {
				return false;
			}
		}
	}

	return true;
}

/// \brief Basic insertion operator to output a rough summary of an alignment to an ostream
///
/// \relates alignment
ostream & cath::align::operator<<(ostream         &arg_ostream,  ///< The stream to which to output
                                  const alignment &arg_alignment ///< The alignment to summarise
                                  ) {
	const alignment::size_type length      = arg_alignment.length();
	const alignment::size_type num_entries = arg_alignment.num_entries();

	arg_ostream << "alignment[";
	arg_ostream << length << " positions: ";
	for (alignment::size_type index_ctr = 0; index_ctr < length; ++index_ctr) {
		arg_ostream << ((index_ctr > 0) ? "; " : "");
		for (alignment::size_type entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
			const aln_posn_opt position = arg_alignment.position_of_entry_of_index( entry_ctr, index_ctr );
			arg_ostream << ( ( entry_ctr > 0 ) ? " <-> " : "" );
			arg_ostream << ( position ? lexical_cast<string>( *position ) : "" );
		}
	}
	arg_ostream << "]";
	return arg_ostream;
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_vec cath::align::entries_present_at_index(const alignment            &arg_alignment, ///< TODOCUMENT
                                               const alignment::size_type &arg_index      ///< TODOCUMENT
                                               ) {
	return copy_build<size_vec>(
		irange( 0_z, arg_alignment.num_entries() )
			| filtered(
				[&] (const size_t &x) {
					return has_position_of_entry_of_index( arg_alignment, x, arg_index );
				}
			)
	);
}

/// \brief Return the entries that have present positions within the specified range in the specified alignment
///
/// \relates alignment
size_vec cath::align::entries_present_in_index_range(const alignment &arg_alignment,   ///< The alignment to be inspected
                                                     const size_t    &arg_begin_index, ///< Begin index of the range to be inspected
                                                     const size_t    &arg_end_index    ///< One-past-end index of the range to be inspected
                                                     ) {
	return copy_build<size_vec>(
		irange( 0_z, arg_alignment.num_entries() )
			| filtered(
				[&] (const size_t &x) {
					return has_positions_of_entry_in_index_range( arg_alignment, x, arg_begin_index, arg_end_index );
				}
			)
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_t cath::align::num_present_positions_of_index(const alignment            &arg_alignment, ///< TODOCUMENT
                                                   const alignment::size_type &arg_index      ///< TODOCUMENT
                                                   ) {
	return entries_present_at_index( arg_alignment, arg_index ).size();
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_vec cath::align::num_present_positions_by_index(const alignment &arg_alignment ///< TODOCUMENT
                                                     ) {
	const size_t length = arg_alignment.length();
	size_vec num_present_positions;
	num_present_positions.reserve( length );
	for (size_t index_ctr = 0; index_ctr < length; ++index_ctr) {
		num_present_positions.push_back( num_present_positions_of_index( arg_alignment, index_ctr ) );
	}
	return num_present_positions;
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_vec cath::align::present_positions_of_entry(const alignment            &arg_alignment, ///< TODOCUMENT
                                                 const alignment::size_type &arg_entry      ///< TODOCUMENT
                                                 ) {
	const alignment::size_type aln_length = arg_alignment.length();
	size_vec present_positions;
	for ( alignment::size_type index = 0; index < aln_length; ++index ) {
		if ( has_position_of_entry_of_index( arg_alignment, arg_entry, index ) ) {
			present_positions.push_back( index );
		}
	}
	return present_positions;
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_t cath::align::num_present_positions_of_entry(const alignment            &arg_alignment, ///< TODOCUMENT
                                      const alignment::size_type &arg_entry      ///< TODOCUMENT
                                      ) {
	return present_positions_of_entry( arg_alignment, arg_entry ).size();
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_vec cath::align::num_present_positions_by_entry(const alignment &arg_alignment ///< TODOCUMENT
                                                     ) {
	const size_t num_entries = arg_alignment.num_entries();
	size_vec num_present_positions;
	num_present_positions.reserve( num_entries );
	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		num_present_positions.push_back( num_present_positions_of_entry( arg_alignment, entry_ctr ) );
	}
	return num_present_positions;
}

/// \brief Convenience function for finding the index of the first present position of the specified entry in the specified alignment
///
/// \returns An optional containing the index of the first present position if there is one or nothing otherwise
///
/// \relates alignment
aln_size_opt cath::align::get_index_of_first_present_position_of_entry(const alignment            &arg_alignment, ///< The alignment to be queried
                                                                       const alignment::size_type &arg_entry      ///< The entry whose first present position is to be found
                                                                       ) {
	const alignment::size_type aln_length = arg_alignment.length();
	for ( alignment::size_type index = 0; index < aln_length; ++index ) {
		if ( has_position_of_entry_of_index( arg_alignment, arg_entry, index ) ) {
			return index;
		}
	}
	return none;
}

/// \brief Convenience function for finding the index of the last present position of the specified entry in the specified alignment
///
/// \returns An optional containing the index of the last present position if there is one or nothing otherwise
///
/// \relates alignment
aln_size_opt cath::align::get_index_of_last_present_position_of_entry(const alignment            &arg_alignment, ///< The alignment to be queried
                                                                      const alignment::size_type &arg_entry      ///< The entry whose last present position is to be found
                                                                      ) {
	const alignment::size_type aln_length = arg_alignment.length();
	for ( alignment::size_type index_offset_1 = aln_length; index_offset_1 > 0; --index_offset_1 ) {
		const alignment::size_type index = index_offset_1 - 1;
		if ( has_position_of_entry_of_index( arg_alignment, arg_entry, index ) ) {
			return index;
		}
	}
	return none;
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_t cath::align::num_present_positions_of_both_entries(const alignment            &arg_alignment, ///< TODOCUMENT
                                                          const alignment::size_type &arg_entry_a,   ///< TODOCUMENT
                                                          const alignment::size_type &arg_entry_b    ///< TODOCUMENT
                                                          ) {
	return numeric_cast<size_t>( count_if(
		irange( 0_z, arg_alignment.length() ),
		[&] (const size_t &x) {
			return has_position_of_both_entries_of_index( arg_alignment, arg_entry_a, arg_entry_b, x );
		}
	) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_t cath::align::num_present_positions_of_all_entries(const alignment &arg_alignment ///< TODOCUMENT
                                                         ) {
	return numeric_cast<size_t>( count_if(
		irange( 0_z, arg_alignment.length() ),
		[&] (const size_t &x) {
			return has_position_of_all_entries_of_index( arg_alignment, x );
		}
	) );
}

/// \brief Convenience function for finding the index of the last present position of both of the specified two entries in the specified alignment
///
/// \returns An optional containing the index of the last present position if there is one or nothing otherwise
///
/// \relates alignment
aln_size_opt cath::align::get_index_of_first_present_position_of_both_entries(const alignment            &arg_alignment, ///< The alignment to be queried
                                                                              const alignment::size_type &arg_entry_a,   ///< The first entry whose last present position is to be found
                                                                              const alignment::size_type &arg_entry_b    ///< The second entry whose last present position is to be found
                                                                              ) {
	const alignment::size_type aln_length = arg_alignment.length();
	for ( alignment::size_type index = 0; index < aln_length; ++index ) {
		if ( has_position_of_both_entries_of_index( arg_alignment, arg_entry_a, arg_entry_b, index ) ) {
			return index;
		}
	}
	return none;
}

/// \brief Convenience function for finding the index of the last present position of both of the specified two entries in the specified alignment
///
/// \returns An optional containing the index of the last present position if there is one or nothing otherwise
///
/// \relates alignment
aln_size_opt cath::align::get_index_of_last_present_position_of_both_entries(const alignment            &arg_alignment, ///< The alignment to be queried
                                                                             const alignment::size_type &arg_entry_a,   ///< The first entry whose last present position is to be found
                                                                             const alignment::size_type &arg_entry_b    ///< The second entry whose last present position is to be found
                                                                             ) {
	const alignment::size_type aln_length = arg_alignment.length();
	for ( alignment::size_type index_offset_1 = aln_length; index_offset_1 > 0; --index_offset_1 ) {
		const alignment::size_type index = index_offset_1 - 1;
		if ( has_position_of_both_entries_of_index( arg_alignment, arg_entry_a, arg_entry_b, index ) ) {
			return index;
		}
	}
	return none;
}

/// \brief Convenience function for finding first present position of the specified entry in the specified alignment
///
/// \returns An optional containing the position of the first present position if there is one or nothing otherwise
///
/// \relates alignment
aln_posn_opt cath::align::get_first_present_position_of_entry(const alignment            &arg_alignment, ///< The alignment to be queried
                                                              const alignment::size_type &arg_entry      ///< The entry whose first present position is to be found
                                                              ) {
	const aln_size_opt index_of_first_present = get_index_of_first_present_position_of_entry(
		arg_alignment,
		arg_entry
	);
	return ( index_of_first_present ) ? get_position_of_entry_of_index( arg_alignment, arg_entry, index_of_first_present.get() )
	                                  : aln_posn_opt();
}

/// \brief Convenience function for finding last present position of the specified entry in the specified alignment
///
/// \returns An optional containing the position of the last present position if there is one or nothing otherwise
///
/// \relates alignment
aln_posn_opt cath::align::get_last_present_position_of_entry(const alignment            &arg_alignment, ///< The alignment to query
                                                             const alignment::size_type &arg_entry      ///< The entry whose last present position is to be found
                                                             ) {
	const aln_size_opt index_of_last_present = get_index_of_last_present_position_of_entry(
		arg_alignment,
		arg_entry
	);
	return ( index_of_last_present ) ? get_position_of_entry_of_index( arg_alignment, arg_entry, index_of_last_present.get() )
	                                 : aln_posn_opt();
}

/// \brief TODOCUMENT
///
/// \relates alignment
aln_posn_opt cath::align::get_max_last_present_position(const alignment &arg_alignment ///< TODOCUMENT
                                                        ) {
	// Grab the number of entries and sanity check that it isn't zero
	const size_t num_entries = arg_alignment.num_entries();
	if ( num_entries == 0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate get_max_last_present_position() of alignment with no entries"));
	}

	 // Loop over the entries to find the max_last_present_position
	 aln_posn_opt max_last_present_position;
	 for (const auto &entry_ctr : irange( 0_z, num_entries ) ) {
	 	max_last_present_position = max(
	 		max_last_present_position,
	 		get_last_present_position_of_entry( arg_alignment, entry_ctr )
	 	);
	 }
	 return max_last_present_position;
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_type cath::align::get_score_of_entry_and_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                           const size_t    &arg_entry,     ///< TODOCUMENT
                                                           const size_t    &arg_index      ///< TODOCUMENT
                                                           ) {
	const alignment_residue_scores &scores = arg_alignment.get_alignment_residue_scores();
	return get_score_to_present_entries_of_index( scores, arg_entry, arg_index );
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_type cath::align::get_mean_score_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                      const size_t    &arg_index      ///< TODOCUMENT
                                                      ) {
	const size_t num_present_entries = num_present_positions_of_index( arg_alignment, arg_index );
	return get_total_score_of_index( arg_alignment, arg_index ) / numeric_cast<float_score_type>( num_present_entries );
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_type cath::align::get_total_score_of_index(const alignment &arg_alignment, ///< TODOCUMENT
                                                       const size_t    &arg_index      ///< TODOCUMENT
                                                       ) {
	const alignment_residue_scores &scores = arg_alignment.get_alignment_residue_scores();
	const size_vec  present_entries        = entries_present_at_index( arg_alignment, arg_index );
	float_score_vec scores_of_index;
	scores_of_index.reserve( present_entries.size() );
	for (const size_t &present_entry : present_entries) {
		scores_of_index.push_back( get_score_to_present_entries_of_index( scores, present_entry, arg_index ) );
	}

	return accumulate( scores_of_index, 0.0 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_vec cath::align::get_mean_score_by_index(const alignment &arg_alignment ///< TODOCUMENT
                                                     ) {
	const size_t length = arg_alignment.length();
	float_score_vec mean_score_by_index;
	mean_score_by_index.reserve( length );
	for (size_t index_ctr = 0; index_ctr < length; ++index_ctr) {
		mean_score_by_index.push_back( get_mean_score_of_index( arg_alignment, index_ctr ) );
	}
	return mean_score_by_index;
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_vec cath::align::get_total_score_by_index(const alignment &arg_alignment ///< TODOCUMENT
                                                      ) {
	const size_t length = arg_alignment.length();
	float_score_vec total_score_by_index;
	total_score_by_index.reserve( length );
	for (size_t index_ctr = 0; index_ctr < length; ++index_ctr) {
		total_score_by_index.push_back( get_total_score_of_index( arg_alignment, index_ctr ) );
	}
	return total_score_by_index;
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_vec cath::align::get_total_score_or_num_positions_by_index(const alignment &arg_alignment ///< TODOCUMENT
                                                                       ) {
	if ( arg_alignment.is_scored() ) {
		return get_total_score_by_index( arg_alignment );
	}
	const size_vec num_present_positions = num_present_positions_by_index(arg_alignment);
	float_score_vec num_present_positions_float;
	num_present_positions_float.reserve( num_present_positions.size() );
	for (const size_t &num : num_present_positions) {
		num_present_positions_float.push_back( numeric_cast<float_score_type>( num ) );
	}
	return num_present_positions_float;

}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::set_scores(alignment               &arg_alignment, ///< TODOCUMENT
                             const score_opt_vec_vec &arg_scores     ///< TODOCUMENT
                             ) {
	arg_alignment.set_scores( make_alignment_residue_scores( arg_alignment, arg_scores ) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
alignment cath::align::set_scores_copy(alignment                arg_alignment, ///< TODOCUMENT
                                       const score_opt_vec_vec &arg_scores     ///< TODOCUMENT
									   ) {
	set_scores( arg_alignment, arg_scores );
	return arg_alignment;
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::set_empty_scores(alignment &arg_alignment ///< TODOCUMENT
                                   ) {
	const alignment::size_type length      = arg_alignment.length();
	const alignment::size_type num_entries = arg_alignment.num_entries();
	set_scores(
		arg_alignment,
		score_opt_vec_vec( num_entries, score_opt_vec( length ) )
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
alignment cath::align::set_empty_scores_copy(alignment arg_alignment ///< TODOCUMENT
                                             ) {
	set_empty_scores( arg_alignment);
	return arg_alignment;
}

/// \brief TODOCUMENT
///
/// \relates alignment
alignment cath::align::make_single_alignment(const size_t &arg_length
                                             ) {
	alignment new_alignment( 1 );
	for (size_t index_ctr = 0; index_ctr < arg_length; ++index_ctr) {
		new_alignment.set_position_value( 0, index_ctr, index_ctr );
	}
	return new_alignment;
}

/// \brief TODOCUMENT
///
/// \relates alignment
alignment cath::align::alignment_offset_1_factory(aln_posn_opt_vec_vec arg_data ///< TODOCUMENT
                                                  ) {
	// Remove the offset_1 from (the copy of) the data
	for (aln_posn_opt_vec &data_part : arg_data) {
		for (aln_posn_opt &data_element : data_part) {
			if ( data_element ) {
				size_t &value_ref = *data_element;
				check_offset_1( value_ref );
				--value_ref;
			}
		}
	}
	// Return an alignment
	return alignment( arg_data );
}

/// \brief Build a new alignment by applying the specified permutation to the entries of the specified alignment
///
/// \relates alignment
///
/// \pre The size of the permutation must match the number of entries in the alignment
///      else an invalid_argument_exception will be thrown
///
/// \pre The permutation must be valid permutation (ie contain one copy of
///      each number 0, 1, ..., (num_entries - 1) but not necessarily in that order )
///      else an invalid_argument_exception will be thrown
///
/// \returns An alignment
alignment cath::align::make_permuted_alignment(const alignment &arg_alignment,  ///< The alignment to use a source
                                               const size_vec  &arg_permutation ///< The permutation to apply to the specified alignment
                                               ) {
	// Sanity check the number of entries matches the size of the permutation
	const size_t num_entries      = arg_alignment.num_entries();
	const size_t aln_length       = arg_alignment.length();
	const size_t permutation_size = arg_permutation.size();
	if ( num_entries != permutation_size ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Alignment with "
			+ lexical_cast<string>( num_entries      )
			+ " entries cannot be permuted by permutation of size "
			+ lexical_cast<string>( permutation_size )
		));
	}

	// Construct a new permuted alignment
	alignment new_alignment( num_entries );
	for (size_t index = 0; index < aln_length; ++index) {
		for (size_t entry = 0; entry < num_entries; ++entry) {
			const size_t permuted_entry = arg_permutation[ entry ];
			const aln_posn_opt position = arg_alignment.position_of_entry_of_index( entry, index );
			if ( position ) {
				new_alignment.set_position_value( permuted_entry, index, *position );
			}
		}
	}

	// Return the new alignment
	return new_alignment;
}

/// \brief TODOCUMENT
///
/// \relates alignment
/// \relatesalso alignment_row
void cath::align::append_row(alignment           &arg_alignment, ///< TODOCUMENT
                             const alignment_row &arg_aln_row    ///< TODOCUMENT
                             ) {
	// Sanity check the input
	const alignment::size_type aln_num_entries = arg_alignment.num_entries();
	const alignment::size_type aln_length      = arg_alignment.length();
	const size_t               num_row_entries = arg_aln_row.num_entries();
	if ( aln_num_entries != num_row_entries ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"arg_aln_row with "
			+ lexical_cast<string>( num_row_entries )
			+ " entries cannot be appended to alignment with "
			+ lexical_cast<string>( aln_num_entries )
			+ " entries"
			));
	}

	// Add each of the entries
	for (alignment::size_type entry_ctr = 0; entry_ctr < aln_num_entries; ++entry_ctr) {
		const aln_posn_opt &position = arg_aln_row.position_of_entry( entry_ctr );
		if ( position ) {
			set_position_value( arg_alignment, entry_ctr, aln_length, *position );
		}
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::append_row(alignment              &arg_alignment,          ///< TODOCUMENT
                             const aln_posn_opt_vec &arg_positions           ///< TODOCUMENT
                             ) {
	append_row( arg_alignment, alignment_row( arg_positions) );
}

/// \brief TODOCUMENT
///
/// \relates alignment
void cath::align::set_position_value(alignment                  &arg_alignment, ///< TODOCUMENT
                                     const alignment::size_type &arg_entry,     ///< TODOCUMENT
                                     const alignment::size_type &arg_index,     ///< TODOCUMENT
                                     const size_t               &arg_value      ///< TODOCUMENT
                                     ) {
	arg_alignment.set_position_value(
		arg_entry,
		arg_index,
		arg_value
	);
}

