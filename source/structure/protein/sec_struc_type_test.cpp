/// \file
/// \brief The sec_struc_type test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "exception/invalid_argument_exception.h"
#include "structure/protein/sec_struc_type.h"

#include <sstream>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;

namespace cath {
	namespace test {

		/// \brief The sec_struc_type_test_suite_fixture to assist in testing sec_struc_type
		struct sec_struc_type_test_suite_fixture {
		protected:
			~sec_struc_type_test_suite_fixture() noexcept = default;
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(sec_struc_type_test_suite, cath::test::sec_struc_type_test_suite_fixture)

/// \brief Test that parsing from an "H" works
BOOST_AUTO_TEST_CASE(input_helix) {
	BOOST_CHECK_EQUAL(sec_struc_type::ALPHA_HELIX, lexical_cast<sec_struc_type>("H"));
}

/// \brief Test that parsing from an "S" works
BOOST_AUTO_TEST_CASE(input_strand) {
	BOOST_CHECK_EQUAL(sec_struc_type::BETA_STRAND, lexical_cast<sec_struc_type>("S"));
}

/// \brief Test that parsing from an "nonsense" throws
BOOST_AUTO_TEST_CASE(input_nonsense) {
	BOOST_CHECK_THROW(lexical_cast<sec_struc_type>("nonsense"), invalid_argument_exception);
}

/// \brief Test that outputting an sec_struc_type::ALPHA_HELIX works
BOOST_AUTO_TEST_CASE(output_helix) {
	BOOST_CHECK_EQUAL("H", lexical_cast<string>(sec_struc_type::ALPHA_HELIX));
}

/// \brief Test that outputting an sec_struc_type::BETA_STRAND works
BOOST_AUTO_TEST_CASE(output_strand) {
	BOOST_CHECK_EQUAL("S", lexical_cast<string>(sec_struc_type::BETA_STRAND));
}

/// \brief Test that outputting an sec_struc_type::COIL works
BOOST_AUTO_TEST_CASE(output_coil) {
	BOOST_CHECK_EQUAL(" ", lexical_cast<string>(sec_struc_type::COIL));
}

BOOST_AUTO_TEST_SUITE_END()

