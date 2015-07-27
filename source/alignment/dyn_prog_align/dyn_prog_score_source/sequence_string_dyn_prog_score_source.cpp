/// \file
/// \brief The sequence_string_dyn_prog_score_source class definitions

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

#include "sequence_string_dyn_prog_score_source.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "exception/invalid_argument_exception.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::algorithm::is_upper;
using boost::numeric_cast;

/// \brief Check that the specified string is valid (current checks: all characters are upper-case letters)
///        and throw invalid_argument_exception if not
void sequence_string_dyn_prog_score_source::check_sequence_string(const string &arg_sequence_string ///< The string to be checked
                                                                  ) {
	if ( !all( arg_sequence_string, is_upper() ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Sequence string contains characters other than upper-case letters"));
	}
}

/// \brief Return the number of elements in the first entry to by aligned with dynamic-programming
size_t sequence_string_dyn_prog_score_source::do_get_length_a() const {
	return sequence_string_a.length();
}

/// \brief Return the number of elements in the second entry to by aligned with dynamic-programming
size_t sequence_string_dyn_prog_score_source::do_get_length_b() const {
	return sequence_string_b.length();
}

/// \brief Return 1 if the letters at the specified indices match and are not Xs, or return 0 otherwise
score_type sequence_string_dyn_prog_score_source::do_get_score(const size_t &arg_index_a, ///< The index of the element of interest in the first  sequence
                                                               const size_t &arg_index_b  ///< The index of the element of interest in the second sequence
                                                               ) const {
	const char &char_a = sequence_string_a[arg_index_a];
	const char &char_b = sequence_string_b[arg_index_b];

	const bool chars_match          = (char_a == char_b);
	const bool chars_are_not_both_x = (char_a != 'X' || char_b != 'X');

	return (chars_match && chars_are_not_both_x) ? 1 : 0;
}

/// \brief Ctor for sequence_string_dyn_prog_score_source
sequence_string_dyn_prog_score_source::sequence_string_dyn_prog_score_source(const string &arg_sequence_string_a, ///< The first  sequence to align
                                                                             const string &arg_sequence_string_b  ///< The second sequence to align
                                                                             ) : sequence_string_a ( arg_sequence_string_a ),
                                                                                 sequence_string_b ( arg_sequence_string_b ) {
	// Check that both sequences strings are valid
	check_sequence_string(sequence_string_a);
	check_sequence_string(sequence_string_b);
}

