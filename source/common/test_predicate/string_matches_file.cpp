/// \file
/// \brief The string_matches_file class definitions

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

#include "string_matches_file.hpp"

#include "common/file/open_fstream.hpp"
#include "common/file/slurp.hpp"
#include "common/file/spew.hpp"
#include "common/test_predicate/detail/strings_equal.hpp"
#include "common/test_predicate/files_equal.hpp"

#include <fstream>

using namespace cath;
using namespace cath::common;

using boost::filesystem::path;
using boost::test_tools::predicate_result;
using std::ifstream;
using std::string;

/// \brief Ctor for string_matches_file
string_matches_file::string_matches_file(const bool          &arg_overwrite_diff_expected_with_got, ///< TODOCUMENT
                                         const str_size_type &arg_diff_half_width                   ///< TODOCUMENT
                                         ) : overwrite_diff_expected_with_got { arg_overwrite_diff_expected_with_got },
                                             diff_half_width                  { arg_diff_half_width                  } {
}

/// \brief Ctor for string_matches_file
string_matches_file::string_matches_file(const str_size_type &arg_diff_half_width ///< TODOCUMENT
                                         ) : diff_half_width { arg_diff_half_width } {
}

/// \brief TODOCUMENT
predicate_result string_matches_file::operator()(const string &arg_got_string, ///< TODOCUMENT
                                                 const path   &arg_filename    ///< TODOCUMENT
                                                 ) const {
	const string expected_file = slurp( arg_filename );

	const string got_name{ "got string" };

	const predicate_result the_result = test::detail::strings_equal(
		arg_got_string,
		got_name,
		expected_file,
		files_equal::FILENAME_NAME_PREFIX + arg_filename.string(),
		diff_half_width
	);

	// If overwrite_diff_expected_with_got and the result's negative,
	// overwrite the expected (arg_filename) with the got (arg_got_string)
	if ( ! the_result && ( overwrite_diff_expected_with_got || ( getenv( "BOOTSTRAP_TESTS" ) != nullptr ) ) ) {
		spew( arg_filename, arg_got_string );
	}

	return the_result;
}

