/// \file
/// \brief The files_equal class definitions

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

#include "files_equal.hpp"

#include "common/file/open_fstream.hpp"
#include "test/predicate/istream_and_file_equal.hpp"

#include <fstream>

using namespace cath::common;
using namespace cath::test;

using boost::filesystem::path;
using boost::test_tools::predicate_result;
using std::getenv;
using std::ifstream;
using std::string;

/// \brief TODOCUMENT
const string files_equal::FILENAME_NAME_PREFIX{ "file " };

/// \brief Ctor for files_equal
files_equal::files_equal(const bootstrap_mode &arg_bootstrapping,  ///< TODOCUMENT
                         const str_size_type  &arg_diff_half_width ///< TODOCUMENT
                         ) : bootstrapping   { arg_bootstrapping   },
                             diff_half_width { arg_diff_half_width } {
}

/// \brief Ctor for files_equal
files_equal::files_equal(const str_size_type &arg_diff_half_width ///< TODOCUMENT
                         ) : diff_half_width { arg_diff_half_width } {
}

/// \brief TODOCUMENT
predicate_result files_equal::operator()(const path &arg_filename1, ///< TODOCUMENT
                                         const path &arg_filename2  ///< TODOCUMENT
                                         ) const {
	// Return false if the two files are literally the same
	// because there is no good reason
	if ( arg_filename1 == arg_filename2 ) {
		predicate_result result( false );
		result.message() << "Files are literally both the same file (" << arg_filename1 << ")";
		return result;
	}

	// Otherwise, create an ifstream for the first file...
	ifstream file_ifstream1;
	open_ifstream( file_ifstream1, arg_filename1 );

	// ...and then just use a istream_and_file_equal
	const predicate_result the_result = istream_and_file_equal{ bootstrapping, diff_half_width }(
		file_ifstream1,
		FILENAME_NAME_PREFIX + arg_filename1.string(),
		arg_filename2
	);

	return the_result;
}
