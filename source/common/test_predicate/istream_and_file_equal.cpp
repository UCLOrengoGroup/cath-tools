/// \file
/// \brief The istream_and_file_equal class definitions

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

#include "istream_and_file_equal.hpp"

#include "common/file/open_fstream.hpp"
#include "common/test_predicate/files_equal.hpp"

#include <fstream>

using namespace boost::test_tools;
using namespace cath;
using namespace cath::common;
using namespace std;

using boost::filesystem::path;

/// \brief Ctor for istream_and_file_equal
istream_and_file_equal::istream_and_file_equal(const bool          &arg_overwrite_diff_expected_with_got, ///< TODOCUMENT
                                               const str_size_type &arg_diff_half_width                   ///< TODOCUMENT
                                               ) : overwrite_diff_expected_with_got( arg_overwrite_diff_expected_with_got ),
                                                   diff_half_width(arg_diff_half_width) {
}

/// \brief Ctor for istream_and_file_equal
istream_and_file_equal::istream_and_file_equal(const str_size_type &arg_diff_half_width ///< TODOCUMENT
                                               ) : overwrite_diff_expected_with_got( false               ),
                                                   diff_half_width                 ( arg_diff_half_width ) {
}

/// \brief TODOCUMENT
predicate_result istream_and_file_equal::operator()(istream      &arg_istream, ///< TODOCUMENT
                                                    const string &arg_name,    ///< TODOCUMENT
                                                    const path   &arg_filename ///< TODOCUMENT
                                                    ) const {
	ifstream file_ifstream;
	open_ifstream(file_ifstream, arg_filename);

	istreams_equal the_istreams_equal(diff_half_width);
	const predicate_result the_result = the_istreams_equal(
		arg_istream,
		arg_name,
		file_ifstream,
		files_equal::FILENAME_NAME_PREFIX + arg_filename.string()
	);

	// If overwrite_diff_expected_with_got and the result's negative,
	// overwrite the expected (arg_filename) with the got (arg_istream)
	if ( overwrite_diff_expected_with_got && ! the_result ) {
		ofstream overwrite_ofstream;
		open_ofstream( overwrite_ofstream, arg_filename );
		arg_istream.clear();
		arg_istream.seekg( 0 );
		overwrite_ofstream << arg_istream.rdbuf();
		overwrite_ofstream.close();
	}

	return the_result;
}

