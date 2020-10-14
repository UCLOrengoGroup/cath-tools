/// \file
/// \brief The hits_processor test suite

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

#include <boost/test/unit_test.hpp>

#include "cath/resolve_hits/read_and_process_hits/hits_processor/summarise_hits_processor.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/write_html_hits_processor.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/write_json_hits_processor.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/write_results_hits_processor.hpp"

#include <sstream>

namespace cath { namespace test { } }

using namespace ::cath::rslv::detail;
using namespace ::cath::test;

using ::std::ostream;
using ::std::ostringstream;

namespace cath {
	namespace test {

		/// \brief The hits_processor_test_suite_fixture to assist in testing hits_processor
		struct hits_processor_test_suite_fixture {
		protected:
			~hits_processor_test_suite_fixture() noexcept {
				try {
					BOOST_CHECK_EQUAL( test_ss.str(), "" );
				}
				catch (...) {}
			}

			/// \brief A test ostream to which hits_processors can write
			ostringstream test_ss;

			/// \brief A ref_vec of ostreams to which hits_processors can write
			ref_vec<ostream> ostreams{ { test_ss } };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(hits_processor_test_suite, hits_processor_test_suite_fixture)



BOOST_AUTO_TEST_SUITE(requires_strictly_worse_hits_works)

BOOST_AUTO_TEST_CASE(summarise_hits_processor_requires_strictly_worse_hits) {
	BOOST_CHECK(   summarise_hits_processor    ( ostreams ).requires_strictly_worse_hits() );
}

BOOST_AUTO_TEST_CASE(write_html_hits_processor_requires_strictly_worse_hits) {
	BOOST_CHECK(   write_html_hits_processor   ( ostreams ).requires_strictly_worse_hits() );
}

BOOST_AUTO_TEST_CASE(write_json_hits_processor_does_not_require_strictly_worse_hits) {
	BOOST_CHECK( ! write_json_hits_processor   ( ostreams ).requires_strictly_worse_hits() );
}

BOOST_AUTO_TEST_CASE(write_results_hits_processor_does_not_require_strictly_worse_hits) {
	BOOST_CHECK( ! write_results_hits_processor( ostreams ).requires_strictly_worse_hits() );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()
