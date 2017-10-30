/// \file
/// \brief The string_aligner class definitions

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

#include "string_aligner.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/find_format.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/formatter.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/count.hpp>
#include <boost/range/algorithm/count_if.hpp>

#include "alignment/alignment.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.hpp"
#include "alignment/gap/gap_penalty.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/size_t_literal.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/out_of_range_exception.hpp"

#include <algorithm>
#include <string>

using namespace boost::algorithm;
using namespace cath;
using namespace cath::align::detail;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::is_space;
using boost::algorithm::is_upper;
using boost::algorithm::token_compress_on;
using boost::algorithm::trim_copy;
using boost::lexical_cast;
using boost::numeric_cast;
using boost::range::count;
using boost::range::count_if;

/// \brief An NVI pass-through to the derived class's concrete definition of the virtual do_align() method
str_str_score_tpl string_aligner::align(const string      &arg_sequence_string_a, ///< The first  string to be aligned
                                        const string      &arg_sequence_string_b, ///< The second string to be aligned
                                        const gap_penalty &arg_gap_penalty        ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                        ) const {
	// Pass-through to the virtual do_align() method to perform the aligning
	const str_str_pair aligned_strings = do_align(
		arg_sequence_string_a,
		arg_sequence_string_b,
		arg_gap_penalty
	);

	// Perform some checks on the output
	check_aligned_string_pair_is_valid(    aligned_strings.first,  aligned_strings.second );
	check_aligned_string_matches_original( aligned_strings.first,  arg_sequence_string_a  );
	check_aligned_string_matches_original( aligned_strings.second, arg_sequence_string_b  );

	// Return the two aligned strings, along with a freshly computed score
	return make_tuple(
		aligned_strings.first,
		aligned_strings.second,
		get_score_of_aligned_sequence_strings(
			aligned_strings.first,
			aligned_strings.second,
			arg_gap_penalty
		)
	);
}

/// \brief Check that the sequence of (upper-case) letters in an aligned string matches that in the original string
void cath::align::detail::check_aligned_string_matches_original(const string &arg_aligned_string, ///< The aligned version of the original string
                                                                const string &arg_orig            ///< The original string that was to be aligned
                                                                ) {
	const string letters_from_aligned = find_format_all_copy( arg_aligned_string, token_finder( ! is_upper() ), empty_formatter(arg_aligned_string) );
	if ( arg_orig != letters_from_aligned ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Letters in aligned string do not match those in original"));
	}
}

/// \brief Perform sanity checks on a pair of aligned strings
void cath::align::detail::check_aligned_string_pair_is_valid(const string &arg_aligned_string_a, ///< The aligned version of the first string
                                                             const string &arg_aligned_string_b  ///< The aligned version of the second string
                                                             ) {
	// Check the two strings are of equal length
	if ( arg_aligned_string_a.length() != arg_aligned_string_b.length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Supposedly aligned strings \""
			+ arg_aligned_string_a
			+ "\" and \""
			+ arg_aligned_string_b
			+ "\" are not of equal length"
		));
	}

	// Perform further checks on each string
	check_aligned_string_is_valid( arg_aligned_string_a );
	check_aligned_string_is_valid( arg_aligned_string_b );
}


/// \brief Check that the specified aligned string is valid
///
/// \pre The specified must be a valid aligned string:
///       * it must contain nothing but spaces, upper-case letters and '-' characters
///       * spaces must all be grouped at the end
///       * the central, non-space character must begin and end with upper-case letters
///      else an invalid_argument_exception will be thrown
void cath::align::detail::check_aligned_string_is_valid(const string &arg_aligned_string ///< An aligned version of a string
                                                        ) {
	// Create a trimmed copy with spaces removed from the end
	const string trimmed = trim_copy(arg_aligned_string);

	// Check that the first & last non-space characters are upper-case letters
	if ( ! trimmed.empty() && ( ! is_upper()( trimmed.front() ) || !is_upper()( trimmed.back() ) ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Supposedly aligned string \""
			+ trimmed
			+ "\" begins/ends with a non-upper-case-letter character"
		));
	}

	// Check that trimmed portions contain nothing but upper-case letters and '-' characters for gaps
	if ( ! all( trimmed, is_upper() || is_any_of( "-") ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Supposedly aligned string \""
			+ trimmed
			+ "\" contains characters other than upper-case-letters and '-' in the trimmed portion"
		));
	}
}

/// \brief Count the number of open gaps and additional extensions in the string
///
/// For example:
///  * "A-B"       has one open gap and zero extensions
///  * "A---B"     has one open gap and one  extensions
///  * "A-B-C"     has two open gaps and no extensions
///  * "A---B---C" has two open gaps and four extensions
size_size_pair cath::align::detail::get_num_gaps_and_extensions(const string &arg_aligned_string ///< The aligned string (must contain upper-case letters, possibly split by '-' characters and possibly with a head/tail of space characters)
                                                                ) {
	// Check that this is a sensible aligned string
	check_aligned_string_is_valid(arg_aligned_string);

	// Count the number of gaps
	const str_vec gap_split_parts     = split_build<str_vec>( arg_aligned_string, is_any_of( "-" ), token_compress_on );
	const size_t  num_gap_split_parts = gap_split_parts.size();
	const size_t  num_open_gaps       = max( static_cast<size_t>( 1_z ), num_gap_split_parts ) - 1;

	// Count the total number of gap characters and then subtract the number of gap opens
	// to get the number of extensions
	const size_t num_gap_chars      = numeric_cast<size_t>( count( arg_aligned_string, '-' ) );
	const size_t num_gap_extensions = num_gap_chars - num_open_gaps;

	// Return the pair of calculated values
	return make_pair( num_open_gaps, num_gap_extensions );
}

/// \brief Get the score accumulated by the specified pair of strings with the specified gap penalty
///
/// This uses the identity matrix for the positive scores
score_type cath::align::detail::get_score_of_aligned_sequence_strings(const string      &arg_aligned_string_a, ///< The aligned version of the first string
                                                                      const string      &arg_aligned_string_b, ///< The aligned version of the second string
                                                                      const gap_penalty &arg_gap_penalty       ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                                                      ) {
	// Sanity check the input strings
	check_aligned_string_pair_is_valid( arg_aligned_string_a, arg_aligned_string_b );

	// Count the number of matching characters
	const string::size_type str_length = arg_aligned_string_a.length();
	size_t num_matches = 0;
	for (const string::size_type &char_ctr : indices( str_length ) ) {
		if ( arg_aligned_string_a[ char_ctr ] == arg_aligned_string_b[ char_ctr ] ) {
			++num_matches;
		}
	}

	// Count the total number of gap characters
	//
	/// \todo Consider switching all this to use cath::align::gap::gap_open_and_extend_counts_of_alignment()
	const size_size_pair gap_counts_a      = get_num_gaps_and_extensions( arg_aligned_string_a );
	const size_size_pair gap_counts_b      = get_num_gaps_and_extensions( arg_aligned_string_b );
	const size_t         total_gap_opens   = gap_counts_a.first  + gap_counts_b.first;
	const size_t         total_gap_extends = gap_counts_a.second + gap_counts_b.second;

	// Return the score for the matches less the penalties for the gaps
	constexpr score_type score_per_match = 1;

	return   ( numeric_cast<score_type>( num_matches       ) * score_per_match                          )
	       - ( numeric_cast<score_type>( total_gap_opens   ) * arg_gap_penalty.get_open_gap_penalty()   )
	       - ( numeric_cast<score_type>( total_gap_extends ) * arg_gap_penalty.get_extend_gap_penalty() );
}


///// \brief Use a sequence_string_dyn_prog_score_source to align a pair of sequence strings
/////        with dynamic programming and return the resulting alignment and score.
/////
///// \relates sequence_string_dyn_prog_score_source
//score_alignment_pair cath::align::detail::align_sequence_strings(const dyn_prog_aligner &arg_dyn_prog_aligner,  ///< TODOCUMENT
//                                                                 const string           &arg_sequence_string_a, ///< The first  sequence to align
//                                                                 const string           &arg_sequence_string_b, ///< The second sequence to align
//                                                                 const score_type       &arg_gap_penalty       ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
////                                                                 const size_t     &arg_window_width       ///< The width of the window (details? TODOCUMENT)
//                                                                 ) {
//	const gen_dyn_prog_string_aligner the_string_aligner(arg_dyn_prog_aligner);
//	const str_str_score_tpl bob = the_string_aligner.align(arg_sequence_string_a, arg_sequence_string_b, arg_gap_penalty);
////	const sequence_string_dyn_prog_score_source score_source(
////		arg_sequence_string_a,
////		arg_sequence_string_b
////	);
////	return ssap_code_dyn_prog_aligner().align(
////		score_source,
////		arg_gap_penalty,
////		arg_window_width
////	);
//}

/// \brief Format an alignment of an arbitrary number of sequence strings
///
/// \pre The alignment must contain the same number of entries as arg_sequence_strings.
///
/// \pre The strings must be at least as long as suggested by the maximum positions in the alignment.
///
/// \returns A vector of strings matching the originals but with '-' characters according to the gaps
///          in the alignment.
///
/// \todo Consider moving this elsewhere
///
/// \relates sequence_string_dyn_prog_score_source
str_vec cath::align::detail::format_alignment_strings(const alignment      &arg_alignment,       ///< The alignment to format
                                                      const str_vec &arg_sequence_strings ///< The strings corresponding to the alignment
                                                      ) {
	// Grab the dimensions of the inputs
	const alignment::size_type num_entries          = arg_alignment.num_entries();
	const alignment::size_type num_sequence_strings = arg_sequence_strings.size();
	const alignment::size_type alignment_length     = arg_alignment.length();

	// Check the number of entries in the alignment matches the number of sequence strings
	if ( num_entries != num_sequence_strings ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Number of entries in alignment ("
			+ lexical_cast<string>( num_entries )
			+ ") doesn't match the number of sequence strings ("
			+ lexical_cast<string>( num_sequence_strings )
			+ ")"
		));
	}

	// Prepare a vector of output strings that all contain '-' characters for now
	str_vec formatted_strings( num_entries, string( alignment_length, ' ' ) );

	// Loop over the alignment's length
	for (const size_t &alignment_entry : indices( num_entries ) ) {
		// Grab references to this entry's (input) sequence string and (output) formatted string
		// and grab the input string's length
		string                  &formatted_string  = formatted_strings[ alignment_entry ];
		const string            &sequence_string   = arg_sequence_strings[ alignment_entry ];
		const string::size_type  seq_string_length = sequence_string.length();

		// Loop over each index in the alignment
		for (const size_t &alignment_idx : indices( alignment_length ) ) {
			// Note whether the non-space part of the string has started and whether it's finished
			const bool   has_non_space        =  ! all( formatted_string, is_space() );

			// TODOCUMENT
			const size_t num_leters_so_far    = numeric_cast<size_t>( count_if( formatted_string, is_upper() ) );
			const bool   finished_all_letters = ( num_leters_so_far >= seq_string_length );

			// If this entry has a position at this index, grab it, check it against the (input) sequence string's length
			// and then use it to update the appropriate letter of the (output) formatted string
			const aln_posn_opt position = arg_alignment.position_of_entry_of_index( alignment_entry, alignment_idx );
			if ( position ) {
				if ( *position >= seq_string_length ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception("Position retrieved from alignment is out of range of the corresponding sequence string"));
				}
				formatted_string[alignment_idx] = sequence_string[ *position ];
			}
			else if ( has_non_space && ! finished_all_letters ) {
				formatted_string[alignment_idx] = '-';
			}
		}
	}

	// Return the new formatted strings
	return formatted_strings;
}

/// \brief Format an alignment of a pair of sequence strings
///
/// \pre The alignment must contain two entries.
///
/// \pre The strings must be at least as long as suggested by the maximum positions in the alignment.
///
/// \returns A pair of strings matching the originals but with '-' characters according to the gaps
///          in the alignment.
///
/// \todo Consider moving this elsewhere
///
/// \relates sequence_string_dyn_prog_score_source
str_str_pair cath::align::detail::format_alignment_strings(const alignment &arg_alignment,         ///< The alignment to format
                                                           const string    &arg_sequence_string_a, ///< The first  string corresponding to the alignment
                                                           const string    &arg_sequence_string_b  ///< The second string corresponding to the alignment
                                                           ) {
	const str_vec sequence_strings  = { arg_sequence_string_a,
	                                    arg_sequence_string_b };
	const str_vec formatted_strings = format_alignment_strings(
		arg_alignment,
		sequence_strings
	);
	if ( 2 != formatted_strings.size() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("format_alignment_strings() did not return a pair of strings when given a pair of inputs"));
	}
	return make_pair( formatted_strings[ 0 ], formatted_strings[ 1 ] );
}

///// \brief Align a pair of sequence strings and return a pair of formatted equivalents (with '-'s to indicate gaps)
/////
///// This version allows the window width to be explicitly specified.
/////
///// \pre The sequence strings must contain nothing but upper-case characters
/////
///// \relates sequence_string_dyn_prog_score_source
//str_str_pair cath::align::detail::align_and_format_sequence_strings(const string     &arg_sequence_string_a, ///< The first  sequence to align
//                                                                    const string     &arg_sequence_string_b, ///< The second sequence to align
//                                                                    const score_type &arg_gap_penalty,       ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
//                                                                    const size_t     &arg_window_width       ///< The width of the window (details? TODOCUMENT)
//                                                                    ) {
//	const score_alignment_pair score_and_alignment = align_sequence_strings(
//		arg_sequence_string_a,
//		arg_sequence_string_b,
//		arg_gap_penalty,
//		arg_window_width
//	);
//	return format_alignment_strings(
//		score_and_alignment.second,
//		arg_sequence_string_a,
//		arg_sequence_string_b
//	);
//}

///// \brief Align a pair of sequence strings and return a pair of formatted equivalents (with '-'s to indicate gaps)
/////
///// This version uses a default window width that's large enough to mean that there is no windowing.
/////
///// \pre The sequence strings must contain nothing but upper-case characters
/////
///// \relates sequence_string_dyn_prog_score_source
//str_str_pair cath::align::detail::align_and_format_sequence_strings(const string     &arg_sequence_string_a, ///< The first  sequence to align
//                                                                    const string     &arg_sequence_string_b, ///< The second sequence to align
//                                                                    const score_type &arg_gap_penalty        ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
//                                                                    ) {
//	const size_t length_a   = arg_sequence_string_a.length();
//	const size_t length_b   = arg_sequence_string_b.length();
//	const size_t max_length = max(length_a, length_b);
//	return align_and_format_sequence_strings(
//		arg_sequence_string_a,
//		arg_sequence_string_b,
//		arg_gap_penalty,
//		max_length
//	);
//}
