/// \file
/// \brief The files_equal test suite

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

#include <boost/test/unit_test.hpp>

#include "cath/test/global_test_constants.hpp"
#include "cath/test/predicate/files_equal.hpp"

using namespace cath::test;

using boost::filesystem::path;

namespace cath {
	namespace test {

		// \todo Refactor files_equal_test_suite and istreams_equal_test_suite into one
		//       and make it also test istream_and_file_equal

		/// \brief The files_equal_test_suite_fixture to assist in testing files_equal
		struct files_equal_test_suite_fixture : protected global_test_constants {
		protected:
			~files_equal_test_suite_fixture() noexcept = default;

		public:
			const path compare_file        = { TEST_BASIC_FILE_TEST_DATA_DIR() / "compare_file"        };
			const path compare_file_equal  = { TEST_BASIC_FILE_TEST_DATA_DIR() / "compare_file_equal"  };

			const path compare_file_longer = { TEST_BASIC_FILE_TEST_DATA_DIR() / "compare_file_longer" };
			const path compare_file_diff   = { TEST_BASIC_FILE_TEST_DATA_DIR() / "compare_file_diff"   };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(files_equal_test_suite, cath::test::files_equal_test_suite_fixture)

/// \brief Check that equal files do compare as equal
BOOST_AUTO_TEST_CASE(equal_files_compare_equal) {
	BOOST_CHECK_FILES_EQUAL(     compare_file,        compare_file_equal    );
	BOOST_CHECK_FILES_EQUAL(     compare_file_equal,  compare_file          );

	// \todo Add tests that explicitly check the message
	//
	// In the meantime, you can uncomment these to have a quick look at the message
//	BOOST_CHECK_FILES_EQUAL(     compare_file_diff,   compare_file          );
//	BOOST_CHECK_FILES_EQUAL(     compare_file,        compare_file_diff     );
//	BOOST_CHECK_FILES_EQUAL(     compare_file_longer, compare_file          );
//	BOOST_CHECK_FILES_EQUAL(     compare_file,        compare_file_longer   );

	BOOST_CHECK(  files_equal()( compare_file,        compare_file_equal  ) );
	BOOST_CHECK(  files_equal()( compare_file_equal,  compare_file        ) );
}

/// \brief Check that the same file does not compare as equal
///        (so that it will highlight erroneous comparisons of a file against itself)
BOOST_AUTO_TEST_CASE(the_same_file_does_not_comapare_equal) {
	BOOST_CHECK( !files_equal()( compare_file,        compare_file        ) );
	BOOST_CHECK( !files_equal()( compare_file_equal,  compare_file_equal  ) );
	BOOST_CHECK( !files_equal()( compare_file_diff,   compare_file_diff   ) );
	BOOST_CHECK( !files_equal()( compare_file_longer, compare_file_longer ) );
}

/// \brief Check that a file does not compare equal if it is the same but shorter/longer
BOOST_AUTO_TEST_CASE(longer_and_shorter_files_do_not_compare_equal) {
	BOOST_CHECK( !files_equal( bootstrap_mode::NEVER )( compare_file,        compare_file_longer ) );
	BOOST_CHECK( !files_equal( bootstrap_mode::NEVER )( compare_file_longer, compare_file        ) );
}

/// \brief Check that a file does not compare equal to different file
BOOST_AUTO_TEST_CASE(different_files_do_not_compare_equal) {
	BOOST_CHECK( !files_equal( bootstrap_mode::NEVER )( compare_file,        compare_file_diff   ) );
	BOOST_CHECK( !files_equal( bootstrap_mode::NEVER )( compare_file_diff,   compare_file        ) );
}

BOOST_AUTO_TEST_SUITE_END()

