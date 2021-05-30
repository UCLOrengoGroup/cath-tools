/// \file
/// \brief The alignment gap definitions

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

#include "alignment_gap.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/alignment/gap/gap_penalty.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::gap;
using namespace ::cath::common;
using namespace ::std;

using ::boost::irange;
using ::boost::lexical_cast;
using ::boost::numeric_cast;

/// \brief Count the number of gaps in a pairwise alignment
///
///
/// The best way to do this is to count the number of gaps between.
///
///     [...] x x   x   x   x   x   x   x [...]
///     [...] x   x   x   x   x   x   x x [...]
///
/// Counting gaps this way is a bit simplistic. A more objective method is to count the gaps between each of the pairs
///

/// \brief Return the number of gaps in the specified entry of the specified alignment
///
/// \relates alignment
size_t cath::align::gap::detail::get_naive_num_gaps_of_entry(const alignment &prm_alignment, ///< The alignment in which the gaps should be counted
                                                             const size_t    &prm_entry      ///< The entry within the alignment in which the gaps should be counted
                                                             ) {
	// Get the index of the first and last present position of the entry
	// (and sanity check that they are values for both)
	const aln_size_opt index_of_first_present_position = get_index_of_first_present_position_of_entry(
		prm_alignment,
		prm_entry
	);
	const aln_size_opt index_of_last_present_position  = get_index_of_last_present_position_of_entry(
		prm_alignment,
		prm_entry
	);
	if ( ! index_of_first_present_position || ! index_of_last_present_position ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate the number of gaps of an entry in an alignment with no positions present"));
	}

	// Generate a vector representing of unique presence flags
	vector<bool> uniq_presence_list;
	for (const alignment::size_type &aln_ctr : irange( *index_of_first_present_position, *index_of_last_present_position + 1 ) ) {
		const bool present = has_position_of_entry_of_index( prm_alignment, prm_entry, aln_ctr );
		if (uniq_presence_list.empty() || ( uniq_presence_list.back() != present ) ) {
			uniq_presence_list.push_back(present);
		}
	}

	// If the number of unique presences/absences is odd then something has gone wrong
	if (uniq_presence_list.size() % 2 != 1) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"When counting the number of gaps for entry "
			+ lexical_cast<string>( prm_entry )
			+ " in an alignment, the number of unique presences/absences is "
			+ lexical_cast<string>( uniq_presence_list.size() )
			+ " which is even (but which should be odd)"
		));
	}

	// Return the number of absent/present pairs after the first presence
	return ( ( uniq_presence_list.size() - 1 ) / 2 );
}

/// \brief TODOCUMENT
///
/// \relates alignment
size_size_pair cath::align::gap::detail::gap_open_and_extend_counts_of_pair_in_alignment(const alignment &prm_alignment, ///< TODOCUMENT
                                                                                         const size_t    &prm_entry_a,   ///< TODOCUMENT
                                                                                         const size_t    &prm_entry_b    ///< TODOCUMENT
                                                                                         ) {
	// Initialise the counts to 0 and 0
	size_size_pair open_and_extend_counts = make_pair( 0_z, 0_z );

	// Grab the indices in which these two entries first/last both appear together
	const aln_size_opt first_index = get_index_of_first_present_position_of_both_entries( prm_alignment, prm_entry_a, prm_entry_b );
	const aln_size_opt last_index  = get_index_of_last_present_position_of_both_entries ( prm_alignment, prm_entry_a, prm_entry_b );

	// If there are such indices then look for gaps in between
	if ( first_index && last_index ) {

		// It will be necessary to track when gaps have been opened in a/b, so initialise both to false
		bool gap_open_in_a = false;
		bool gap_open_in_b = false;

		// Loop from the first simultaneous appearance to the last
		for (const size_t &index : irange( *first_index, *last_index ) ) {
			// If this has a position in both entries then close both gaps
			if ( has_position_of_both_entries_of_index( prm_alignment, prm_entry_a, prm_entry_b, index ) ) {
				gap_open_in_a = false;
				gap_open_in_b = false;
			}
			// Else if this has a position in b (and hence a gap in a) then...
			else if ( has_position_of_entry_of_index( prm_alignment, prm_entry_b, index ) ) {
				// If a gap is already in a open then increment the extend count
				if ( gap_open_in_a ) {
					open_and_extend_counts.second++;
				}
				// Else this is creating a gap in a so:
				//  * increment the open count and
				//  * set gap_open_in_a to true
				else {
					open_and_extend_counts.first++;
					gap_open_in_a = true;
				}
			}
			// Else if this has a position in a (and hence a gap in b) then...
			else if ( has_position_of_entry_of_index( prm_alignment, prm_entry_a, index ) ) {
				// If a gap is already in b open then increment the extend count
				if ( gap_open_in_b ) {
					open_and_extend_counts.second++;
				}
				// Else this is creating a gap in b so:
				//  * increment the open count and
				//  * set gap_open_in_b to true
				else {
					open_and_extend_counts.first++;
					gap_open_in_b = true;
				}
			}
			// Otherwise there's no position in entry a or b, so just ignore
		}
	}
	// Return the results of the calculations
	return open_and_extend_counts;
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_type cath::align::gap::gap_count_of_alignment(const alignment &prm_alignment ///< The alignment to assess
                                                          ) {
	return gap_open_and_extend_counts_of_alignment( prm_alignment ).first;
}

/// \brief Calculate open and extend gap counts for an alignment
///
/// This is not just a naive count; for that use get_naive_num_gaps()
///
/// This is uses an average of pairwise gap counts to get around some of the problems of naive counts
///
/// Naive Counts Problem #1
/// -----------------------
///
/// A single residue deletion in one sequence represents a small change, so it's only penalised with a single open:
///
///     x x x
///     x - x
///     x x x
///     x x x
///
/// Yet a single insertion, gets penalised with (n - 1) opens:
///
///     x x - x
///     x x x x
///     x x - x
///     x x - x
///
/// Like the deletion, the insertion is a small change which has left the alignment between all (n-1)(n-2)/2
/// other pairs of sequences unaffected.
///
/// Naive Counts Problem #2
/// -----------------------
///
/// The following alignment appears to have 8 gaps:
///
///     [...] x x - - x [...]
///     [...] x - - x x [...]
///     [...] x - x - x [...]
///     [...] x - x - x [...]
///     [...] x - x - x [...]
///
/// ...but after rearranging (whilst completely preserving which residues are aligned with which):
///
///     [...] x x - - x [...]
///     [...] x - x - x [...]
///     [...] x - - x x [...]
///     [...] x - - x x [...]
///     [...] x - - x x [...]
///
/// ...there are now only 6 gaps.
///
/// A Solution
/// ----------
///
/// Both these problems point to a gap count that depends not on the way the full alignment happens to be laid out,
/// but the actual insertions and deletions that interrupt the alignments between each of the pairs of sequences.
///
/// So for each, sequence, one can sum the gap penalties accrued by
///
/// Counting gaps this way is a bit simplistic. A more objective method is to count the gaps between each of the n(n-1)/2 pairs
/// and then averaging (which in the above case gives 14 / 10 = 1.4).
///
/// This can then be multiplied by the number of entries for a sensible overall value (eg 5*1.4 = 7 in the above case).
/// though this is unfortunately not always an integer, eg the following alignment gets an overall gap count of 8/3 = 2.6666...
///
///     [...] x x x [...]
///     [...] x - x [...]
///     [...] x - x [...]
///     [...] x x x [...]
///
/// 3.5 * opening
/// x x - - x
/// x - x - x
/// x - - x x
/// x - - x x
/// x - - x x
///
///
/// x(n - x) / (n - 1)
///
///
/// x x x
/// x - x
/// x - x
/// x x x
///
/// x x x
/// x - x
/// x - x
/// x x x
///
/// x x x
/// x - x
/// x - x
/// x x x
/// ...there are now only 6 gaps.
///
/// Counting gaps this way is a bit simplistic. A more objective method is to count the gaps between each of the n(n-1)/2 pairs
/// and then averaging (which in the above case gives 14 / 10 = 1.4).
///
/// This can then be multiplied by the number of entries for a sensible overall value (eg 5*1.4 = 7 in the above case).
/// though this is unfortunately not always an integer, eg the following alignment gets an overall gap count of 8/3 = 2.6666...
///
///     [...] x x x [...]
///     [...] x - x [...]
///     [...] x - x [...]
///     [...] x x x [...]

///
/// \relates alignment
float_score_float_score_pair cath::align::gap::gap_open_and_extend_counts_of_alignment(const alignment &prm_alignment ///< The alignment to assess
                                                                                       ) {
	check_entry_positions_are_consecutive( prm_alignment );

	const size_t num_entries = prm_alignment.num_entries();

	// Loop over distinct pairs, to accumulate the open/extend counts for each
	size_size_pair total_open_and_extend_counts = make_pair( 0_z, 0_z );
	for (const size_t &entry_a : indices( num_entries ) ) {
		for (const size_t &entry_b : irange( entry_a, num_entries ) ) {
			// Calculate the counts for the pair and add them to the running count
			const size_size_pair pair_couts = detail::gap_open_and_extend_counts_of_pair_in_alignment(
				prm_alignment,
				entry_a,
				entry_b
			);
			total_open_and_extend_counts.first  += pair_couts.first;
			total_open_and_extend_counts.second += pair_couts.second;
		}
	}

	// Return the (casted) result of dividing by ( num_entries - 1 )
	const auto denominator = numeric_cast<float_score_type>( num_entries - 1 );
	return make_pair(
		numeric_cast<float_score_type>( total_open_and_extend_counts.first  ) / denominator,
		numeric_cast<float_score_type>( total_open_and_extend_counts.second ) / denominator
	);
}

/// \brief TODOCUMENT
///
/// \relates alignment
float_score_type cath::align::gap::gap_penalty_value_of_alignment(const alignment   &prm_alignment,  ///< The alignment to assess
                                                                  const gap_penalty &prm_gap_penalty ///< The gap penalty to apply
                                                                  ) {
	const float_score_float_score_pair open_and_extend_counts = gap_open_and_extend_counts_of_alignment( prm_alignment );
	const float_score_type &open_count   = open_and_extend_counts.first;
	const float_score_type &extend_count = open_and_extend_counts.second;
	return (
		  open_count   * numeric_cast<float_score_type>( prm_gap_penalty.get_open_gap_penalty()   )
		+ extend_count * numeric_cast<float_score_type>( prm_gap_penalty.get_extend_gap_penalty() )
	);
}

/// \brief Return the total number of gaps in the specified alignment
///
/// \relates alignment
size_t cath::align::gap::get_naive_num_gaps(const alignment &prm_alignment ///< The alignment in which the gaps should be counted
                                            ) {
	// Loop over the entries in the alignment and add up the number of gaps of each
	size_t num_gaps = 0;
	const alignment::size_type num_entries = prm_alignment.num_entries();
	for (const size_t &entry_ctr : indices( num_entries ) ) {
		num_gaps += detail::get_naive_num_gaps_of_entry( prm_alignment, entry_ctr );
	}

	// Return the total number of gaps
	return num_gaps;
}

