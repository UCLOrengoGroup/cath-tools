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

#include <filesystem>
#include <fstream>

#include "cath/common/file/open_fstream.hpp"
#include "cath/common/file/spew.hpp"
#include "cath/test/predicate/files_equal.hpp"

using namespace ::cath::common;
using namespace ::cath::test;

using ::boost::test_tools::predicate_result;
using ::std::filesystem::path;
using ::std::ifstream;
using ::std::istream;
using ::std::string;

/// \brief Ctor for istream_and_file_equal
istream_and_file_equal::istream_and_file_equal(const bootstrap_mode &prm_bootstrapping,  ///< TODOCUMENT
                                               const str_size_type  &prm_diff_half_width ///< TODOCUMENT
                                               ) : bootstrapping   { prm_bootstrapping   },
                                                   diff_half_width { prm_diff_half_width } {
}

/// \brief Ctor for istream_and_file_equal
istream_and_file_equal::istream_and_file_equal(const str_size_type &prm_diff_half_width ///< TODOCUMENT
                                               ) : diff_half_width{ prm_diff_half_width } {
}

/// \brief TODOCUMENT
predicate_result istream_and_file_equal::operator()(istream      &prm_istream, ///< TODOCUMENT
                                                    const string &prm_name,    ///< TODOCUMENT
                                                    const path   &prm_filename ///< TODOCUMENT
                                                    ) const {
	ifstream file_ifstream = open_ifstream( prm_filename );

	istreams_equal the_istreams_equal(diff_half_width);
	const predicate_result the_result = the_istreams_equal(
		prm_istream,
		prm_name,
		file_ifstream,
		files_equal::FILENAME_NAME_PREFIX + prm_filename.string()
	);

	// If `should_overwrite( bootstrapping )` and the result is negative,
	// overwrite the expected (prm_filename) with the got (prm_istream)
	if ( ! the_result && should_overwrite( bootstrapping ) ) {
		spew( prm_filename, prm_istream );
	}

	return the_result;
}

