/// \file
/// \brief The prc_scores_entry test suite

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

#include "file/prc_scores_file/prc_scores_entry.hpp"
#include "test/test_tools.hpp"

#include <vector>

//using namespace cath::common;
using namespace cath::common::test;
using namespace cath::file;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The prc_scores_entry_test_suite_fixture to assist in testing prc_scores_entry
		struct prc_scores_entry_test_suite_fixture {
		protected:
			~prc_scores_entry_test_suite_fixture() noexcept = default;

			const prc_scores_entry eg_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00", 15, 209, 219, 25.3, 16.3, 1.6e-11 };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(prc_scores_entry_test_suite, cath::test::prc_scores_entry_test_suite_fixture)

BOOST_AUTO_TEST_CASE(equality_works) {
	check_equality_operators_on_diff_vals_range( vector<prc_scores_entry>{
		eg_entry,

		prc_scores_entry{ "belly",   4, 199, 201, 1, "3cazA00", 15, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 0, 199, 201, 1, "3cazA00", 15, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4,   0, 201, 1, "3cazA00", 15, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199,   0, 1, "3cazA00", 15, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 0, "3cazA00", 15, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "ship",    15, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00",  0, 209, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00", 15,   0, 219, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00", 15, 209,   0, 25.3, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00", 15, 209, 219,  0.0, 16.3, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00", 15, 209, 219, 25.3,  0.0, 1.6e-11 },
		prc_scores_entry{ "1i4dA00", 4, 199, 201, 1, "3cazA00", 15, 209, 219, 25.3, 16.3, 0.0     }
	} );
}

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL(
		to_string( eg_entry ),
		"prc_scores_entry[1i4dA00, 4, 199, 201, 1, 3cazA00, 15, 209, 219, 25.300000, 16.300000, 0.000000]"
	);
}

BOOST_AUTO_TEST_CASE(parses_from_line) {
	BOOST_CHECK_EQUAL(
		prc_scores_entry_from_line( "1i4dA00 4       199     201     1       3cazA00 15      209     219       25.3    16.3   1.6e-11" ),
		eg_entry
	);
}

BOOST_AUTO_TEST_CASE(getters) {
	BOOST_CHECK_EQUAL( eg_entry.get_name_1(),   "1i4dA00" );
	BOOST_CHECK_EQUAL( eg_entry.get_start_1(),    4       );
	BOOST_CHECK_EQUAL( eg_entry.get_end_1(),    199       );
	BOOST_CHECK_EQUAL( eg_entry.get_length_1(), 201       );
	BOOST_CHECK_EQUAL( eg_entry.get_hit_num(),    1       );
	BOOST_CHECK_EQUAL( eg_entry.get_name_2(),   "3cazA00" );
	BOOST_CHECK_EQUAL( eg_entry.get_start_2(),   15       );
	BOOST_CHECK_EQUAL( eg_entry.get_end_2(),    209       );
	BOOST_CHECK_EQUAL( eg_entry.get_length_2(), 219       );
	BOOST_CHECK_EQUAL( eg_entry.get_simple(),    25.3     );
	BOOST_CHECK_EQUAL( eg_entry.get_reverse(),   16.3     );
	BOOST_CHECK_EQUAL( eg_entry.get_evalue(),     1.6e-11 );
}

BOOST_AUTO_TEST_SUITE_END()

