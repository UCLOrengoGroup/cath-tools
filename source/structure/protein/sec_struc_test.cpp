/// \file
/// \brief The sec_struc test suite

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

#include <boost/test/auto_unit_test.hpp>

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "common/boost_check_no_throw_diag.h"
#include "common/type_aliases.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "test/global_test_constants.h"

#include <ostream>

namespace cath { namespace common { class invalid_argument_exception; } }

using namespace cath;
using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The sec_struc_test_suite_fixture to assist in testing sec_struc
		struct sec_struc_test_suite_fixture : protected global_test_constants {
		protected:
			~sec_struc_test_suite_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(sec_struc_test_suite, cath::test::sec_struc_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(sec_struc_planar_angles_throws_on_invalid) {
	const doub_vec &invalid_doubles = INVALID_DOUBLES();
	for (const double &invalid_double : invalid_doubles) {
		BOOST_CHECK_THROW(         sec_struc_planar_angles( invalid_double,            0.0,            0.0 ), invalid_argument_exception );
		BOOST_CHECK_THROW(         sec_struc_planar_angles(            0.0, invalid_double,            0.0 ), invalid_argument_exception );
		BOOST_CHECK_THROW(         sec_struc_planar_angles(            0.0,            0.0, invalid_double ), invalid_argument_exception );
		BOOST_CHECK_NO_THROW_DIAG( sec_struc_planar_angles(            0.0,            0.0,            0.0 )                             );
	}
}



BOOST_AUTO_TEST_SUITE_END()
