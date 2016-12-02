/// \file
/// \brief The istreams_equal class definitions

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#include "istreams_equal.hpp"

#include <boost/numeric/conversion/cast.hpp>

#include "exception/invalid_argument_exception.hpp"

#include <algorithm>
#include <iostream> // ***** TEMPORARY *****
#include <vector>

using namespace boost::test_tools;
using namespace cath;
using namespace cath::common;
using namespace std;

constexpr str_size_type istreams_equal::DEFAULT_DIFF_HALF_WIDTH;

/// \brief TODOCUMENT
str_size_type istreams_equal::index_of_first_difference(const string &arg_string1, ///< TODOCUMENT
                                                        const string &arg_string2  ///< TODOCUMENT
                                                        ) {
	// Step through the characters up to the end of the shorter string
	// and return the index of the first difference found, if any
	const str_size_type min_length = min(arg_string1.length(), arg_string2.length());
	for (str_size_type char_ctr = 0; char_ctr < min_length; ++char_ctr) {
		if (arg_string1[char_ctr] != arg_string2[char_ctr]) {
			return char_ctr;
		}
	}

	// No difference has been found up to the length of the shorter string,
	// so if the longer string is the same length, throw an error
	const str_size_type max_length = max(arg_string1.length(), arg_string2.length());
	if (min_length == max_length) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot find any differences in index_of_first_difference"));
	}

	// Otherwise, the end of the shorter string is the first difference,
	// so return its length
	return min_length;
}

/// \brief Ctor for istreams_equal
istreams_equal::istreams_equal(const str_size_type &arg_diff_half_width ///< TODOCUMENT
                               ) : diff_half_width(arg_diff_half_width) {
}

/// \brief TODOCUMENT
predicate_result istreams_equal::operator()(istream       &arg_istream1, ///< TODOCUMENT
                                            const string  &arg_name1,    ///< TODOCUMENT
                                            istream       &arg_istream2, ///< TODOCUMENT
                                            const string  &arg_name2     ///< TODOCUMENT
                                            ) const {
	// Suck the two istreams into strings
	const string input_string1(
		(istreambuf_iterator<char>(arg_istream1)),
		istreambuf_iterator<char>()
	);
	const string input_string2(
		(istreambuf_iterator<char>(arg_istream2)),
		istreambuf_iterator<char>()
	);

	// If the two strings are equal, return a predicate_result of true
	if (input_string1 == input_string2) {
		return predicate_result(true);
	}
	// Otherwise, attempt to output some diagnostic information about the differences
	predicate_result result(false);
	result.message() << arg_name1 << " and " << arg_name2 << " differ:\n";

	const str_size_type first_diff_idx  = index_of_first_difference(input_string1, input_string2);
	const str_size_type start_of_window =   (  first_diff_idx > diff_half_width )
	                                      ? (  first_diff_idx - diff_half_width )
	                                      : 0;
	const vector<const string *> string_ptrs = { &input_string1,
	                                             &input_string2 };
	for (const string *string_ptr : string_ptrs) {
		const str_size_type string_length  = string_ptr->length();
		const str_size_type end_of_window  = min( first_diff_idx + diff_half_width, string_length );
		result.message() << "\"" << (start_of_window > 0 ? "[...]" : "");
		for (str_size_type char_ctr = start_of_window; char_ctr < end_of_window; ++char_ctr) {
			result.message() << string_ptr->operator[](char_ctr);
		}
		result.message() << (end_of_window < string_length ? "[...]" : "") << "\"\n";
	}
	return result;
}
