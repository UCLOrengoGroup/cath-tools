/// \file
/// \brief The roled_scan_stride test suite

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

#include "scan/detail/stride/roled_scan_stride.hpp"

//#include "test/global_test_constants.hpp"

using namespace cath::scan;
using namespace cath::scan::detail;
//using namespace std;

namespace cath {
	namespace test {

		/// \brief The roled_scan_stride_test_suite_fixture to assist in testing roled_scan_stride
		struct roled_scan_stride_test_suite_fixture {
		protected:
			~roled_scan_stride_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(roled_scan_stride_test_suite, cath::test::roled_scan_stride_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ctor_and_getters_work) {
	BOOST_CHECK( true );
	static_assert( get_this_from_stride ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } ) == 1, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_this_to_stride   ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } ) == 2, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_other_from_stride( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } ) == 3, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_other_to_stride  ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } ) == 4, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_this_from_stride ( roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } ) == 3, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_this_to_stride   ( roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } ) == 4, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_other_from_stride( roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } ) == 1, "roled_scan_stride's ctor and getter should return the relevant value back" );
	static_assert( get_other_to_stride  ( roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } ) == 2, "roled_scan_stride's ctor and getter should return the relevant value back" );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(equality_comparison_works) {
	BOOST_CHECK( true );
	static_assert(     roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } },   "scan_stride's operator== should work" );
	static_assert(     roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } == roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } },   "scan_stride's operator== should work" );

	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } == roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::INDEX, scan_stride{ 1, 2, 3, 4 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 3, 4 } } ), "scan_stride's operator== should work" );

	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 2, 1, 1, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 1, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 2, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 2 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 2, 1, 1, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 2, 1, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 2, 1 } } ), "scan_stride's operator== should work" );
	static_assert( ! ( roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 1 } } == roled_scan_stride{ scan_role::QUERY, scan_stride{ 1, 1, 1, 2 } } ), "scan_stride's operator== should work" );
}

BOOST_AUTO_TEST_SUITE_END()
