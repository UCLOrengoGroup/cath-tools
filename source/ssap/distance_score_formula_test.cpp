/// \file
/// \brief The  test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "ssap/distance_score_formula.hpp"

#include "test/global_test_constants.hpp"

using namespace cath;
using namespace std;

using boost::lexical_cast;

namespace cath {
	namespace test {

		/// \brief The distance_score_formula_test_suite_fixture to assist in testing distance_score_formula
		struct distance_score_formula_test_suite_fixture : protected global_test_constants {
		protected:
			~distance_score_formula_test_suite_fixture() noexcept = default;
		};
	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(distance_score_formula_test_suite, cath::test::distance_score_formula_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( name_of_distance_score_formula::get().size(), num_distance_score_formulae );
	for (const distance_score_formula &value : all_distance_score_formulae) {
		BOOST_CHECK_GT( name_of_distance_score_formula::get().at( value ).length(), 0 );
		BOOST_CHECK_GT( lexical_cast<string>( value ).length(), 0 );
	}
}

BOOST_AUTO_TEST_SUITE_END()
