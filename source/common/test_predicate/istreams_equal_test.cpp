/// \file
/// \brief The istreams_equal test suite

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

#include <boost/test/auto_unit_test.hpp>

#include "common/test_predicate/istreams_equal.hpp"

using namespace cath;

using std::istringstream;
using std::string;

namespace cath {
	namespace test {

		/// \brief The istreams_equal_test_suite_fixture to assist in testing istreams_equal
		struct istreams_equal_test_suite_fixture {
		protected:
			~istreams_equal_test_suite_fixture() noexcept = default;
			istreams_equal_test_suite_fixture();
			
			void reset_stringstreams();

			istringstream test_ss;
			istringstream compare_ss;
			istringstream compare_ss_equal;
			istringstream compare_ss_longer;
			istringstream compare_ss_diff;

			const string test_str           = "This is Tony's test string";

			const string compare_str        = "This is a string that will be compared";
			const string compare_str_equal  = "This is a string that will be compared";

			const string compare_str_longer = "This is a string that will be compared and this version has extra stuff at the end";
			const string compare_str_diff   = "This is a string WOW TEXT INSERTED HERE WOW that will be compared";
		};

		/// \brief TODOCUMENT
		istreams_equal_test_suite_fixture::istreams_equal_test_suite_fixture() {
			reset_stringstreams();
		}
		
		/// \brief TODOCUMENT
		void istreams_equal_test_suite_fixture::reset_stringstreams() {
			compare_ss.str        ( compare_str        );
			compare_ss_equal.str  ( compare_str_equal  );
			compare_ss_longer.str ( compare_str_longer );
			compare_ss_diff.str   ( compare_str_diff   );
		}

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(istreams_equal_test_suite, cath::test::istreams_equal_test_suite_fixture)

/// \brief Check that equal istreams do compare as equal
BOOST_AUTO_TEST_CASE(equal_files_compare_equal) {
	BOOST_CHECK_ISTREAMS_EQUAL(     compare_ss,        "istream1", compare_ss_equal,  "istream2"   );
	reset_stringstreams();
	BOOST_CHECK_ISTREAMS_EQUAL(     compare_ss_equal,  "istream1", compare_ss,        "istream2"   );
	reset_stringstreams();
	BOOST_CHECK(  istreams_equal()( compare_ss,        "istream1", compare_ss_equal,  "istream2" ) );
	reset_stringstreams();
	BOOST_CHECK(  istreams_equal()( compare_ss_equal,  "istream1", compare_ss,        "istream2" ) );
}

/// \brief Check that a istream does not compare equal if it is the same but shorter/longer
BOOST_AUTO_TEST_CASE(longer_and_shorter_files_do_not_compare_equal) {
	BOOST_CHECK( !istreams_equal()( compare_ss,        "istream1", compare_ss_longer, "istream2" ) );
	reset_stringstreams();
	BOOST_CHECK( !istreams_equal()( compare_ss_longer, "istream1", compare_ss,        "istream2" ) );
}

/// \brief Check that a istream does not compare equal to different istream
BOOST_AUTO_TEST_CASE(different_files_do_not_compare_equal) {
	BOOST_CHECK( !istreams_equal()( compare_ss,        "istream1", compare_ss_diff,   "istream2" ) );
	reset_stringstreams();
	BOOST_CHECK( !istreams_equal()( compare_ss_diff,   "istream1", compare_ss,        "istream2" ) );
}

BOOST_AUTO_TEST_SUITE_END()

