/// \file
/// \brief The strings_equal class definitions

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

#include "strings_equal.hpp"

#include <boost/range/irange.hpp>

#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"

using namespace cath::common;
using namespace cath::test::detail;

using boost::irange;
using boost::test_tools::predicate_result;
using std::min;
using std::string;

/// \brief TODOCUMENT
size_t cath::test::detail::index_of_first_difference(const string &prm_string1, ///< TODOCUMENT
                                                     const string &prm_string2  ///< TODOCUMENT
                                                     ) {
	// Step through the characters up to the end of the shorter string
	// and return the index of the first difference found, if any
	const size_t min_length = min( prm_string1.length(), prm_string2.length() );
	for (const size_t &char_ctr : indices( min_length ) ) {
		if ( prm_string1[ char_ctr ] != prm_string2[ char_ctr ] ) {
			return char_ctr;
		}
	}

	// No difference has been found up to the length of the shorter string,
	// so if the longer string is the same length, throw an error
	if ( prm_string1.length() == prm_string2.length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot find any differences in index_of_first_difference"));
	}

	// Otherwise, the end of the shorter string is the first difference,
	// so return its length
	return min_length;
}


/// \brief TODOCUMENT
predicate_result cath::test::detail::strings_equal(const string &prm_string1,        ///< TODOCUMENT
                                                   const string &prm_name1,          ///< TODOCUMENT
                                                   const string &prm_string2,        ///< TODOCUMENT
                                                   const string &prm_name2,          ///< TODOCUMENT
                                                   const size_t &prm_diff_half_width ///< TODOCUMENT
                                                   ) {
	// If the two strings are equal, return a predicate_result of true
	if ( prm_string1 == prm_string2 ) {
		return { true };
	}
	// Otherwise, attempt to output some diagnostic information about the differences
	predicate_result result{ false };
	result.message() << prm_name1 << " and " << prm_name2 << " differ:\n";

	const size_t first_diff_idx  = index_of_first_difference( prm_string1, prm_string2 );
	const size_t start_of_window =   (  first_diff_idx > prm_diff_half_width )
	                               ? (  first_diff_idx - prm_diff_half_width )
	                               : 0;
	const auto string_refs = { cref( prm_string1 ),
	                           cref( prm_string2 ) };
	for (const auto &string_ref : string_refs) {
		const size_t string_length  = string_ref.get().length();
		const size_t end_of_window  = min( first_diff_idx + prm_diff_half_width, string_length );
		result.message() << "\"" << (start_of_window > 0 ? "[...]" : "");
		for (const size_t &char_ctr : irange( start_of_window, end_of_window ) ) {
			result.message() << string_ref.get()[ char_ctr ];
		}
		result.message() << ( end_of_window < string_length ? "[...]" : "" ) << "\"\n";
	}
	return result;
}
