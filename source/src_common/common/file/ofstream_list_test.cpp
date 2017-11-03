/// \file
/// \brief The ofstream_list test suite

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

#include "common/file/ofstream_list.hpp"
#include "common/file/read_string_from_file.hpp"
// #include "common/file/simple_file_read_write.hpp"
#include "common/file/temp_file.hpp"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::test;

using std::ostringstream;
using std::string;

namespace cath {
	namespace test {

		/// \brief The ofstream_list_test_suite_fixture to assist in testing ofstream_list
		struct ofstream_list_test_suite_fixture {
		protected:
			~ofstream_list_test_suite_fixture() noexcept = default;

		public:
			const temp_file test_out_1{ ".ofstream_list_test.temp_file.%%%%-%%%%-%%%%-%%%%" };
			const temp_file test_out_2{ ".ofstream_list_test.temp_file.%%%%-%%%%-%%%%-%%%%" };
			ostringstream   test_ostream;
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(ofstream_list_test_suite, ofstream_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	ofstream_list the_list( test_ostream, "-" );

	const ostream_ref_vec result = the_list.open_ofstreams( {
		get_filename( test_out_1 ),
		"-",
		get_filename( test_out_2 )
	} );

	BOOST_REQUIRE_EQUAL( result.size(), 3 );

	result[ 0 ].get() << "one";
	result[ 1 ].get() << "two";
	result[ 2 ].get() << "three";

	// get_filename
	the_list.close_all();

	BOOST_CHECK_EQUAL( read_string_from_file( get_filename( test_out_1 ) ), "one"   );
	BOOST_CHECK_EQUAL( test_ostream.str(),                                  "two"   );
	BOOST_CHECK_EQUAL( read_string_from_file( get_filename( test_out_2 ) ), "three" );
}

BOOST_AUTO_TEST_SUITE_END()

