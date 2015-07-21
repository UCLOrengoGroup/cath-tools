/// \file
/// \brief The plate_scan test suite

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include <boost/test/auto_unit_test.hpp>

#include "structure/view_cache/detail/plate/plate_scan.h"

namespace cath { namespace index { namespace detail { class plate; } } }

using namespace cath::index::detail;

namespace cath {
	namespace test {

		/// \brief The plate_scan_test_suite_fixture to assist in testing plate_scan
		struct plate_scan_test_suite_fixture {
		protected:
			~plate_scan_test_suite_fixture() noexcept = default;
			
		public:
			void test_plate_scan(const size_t &,
			                     const size_t &,
			                     const size_t &,
			                     const size_t &,
			                     const size_t &,
			                     const size_t &,
			                     const size_t &,
			                     const size_t &) const;
		};

	}
}

///// \brief TODOCUMENT
//void plate_scan_test_suite_fixture::test_plate_scan(const size_t &arg_from_index,     ///< TODOCUMENT
//                                                    const size_t &arg_to_index,       ///< TODOCUMENT
//                                                    const size_t &arg_from_size,      ///< TODOCUMENT
//                                                    const size_t &arg_to_size,        ///< TODOCUMENT
//                                                    const size_t &arg_from_step_size, ///< TODOCUMENT
//                                                    const size_t &arg_to_step_size,   ///< TODOCUMENT
//                                                    const size_t &arg_from_offset,    ///< TODOCUMENT
//                                                    const size_t &arg_to_offset       ///< TODOCUMENT
//                                                    ) const {
//	arg_from_index;
//}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(plate_scan_test_suite, cath::test::plate_scan_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( 0, 0 );
}

BOOST_AUTO_TEST_SUITE_END()
