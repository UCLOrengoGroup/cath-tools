/// \file
/// \brief The align_scaffold test suite

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

#include "cath/alignment/io/align_scaffold.hpp"
#include "cath/alignment/test/alignment_fixture.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath::align;
using namespace ::cath::common;

namespace cath {
	namespace test {

		/// \brief The align_scaffold_test_suite_fixture to assist in testing align_scaffold
		struct align_scaffold_test_suite_fixture : protected alignment_fixture {
		protected:
			~align_scaffold_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(align_scaffold_test_suite, cath::test::align_scaffold_test_suite_fixture)

BOOST_AUTO_TEST_CASE(throws_given_zero_scaffold_strings) {
	BOOST_CHECK_THROW( alignment_of_scaffold_lines( { } ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(constructs_aln_a_b_from_scaffold_lines) {
	BOOST_CHECK_EQUAL( alignment_of_scaffold_lines( { "XXXX", "XXX " } ), aln_a_b );
}

BOOST_AUTO_TEST_CASE(constructs_aln_a_b_from_scaffold_lines_with_dashes) {
	BOOST_CHECK_EQUAL( alignment_of_scaffold_lines( { "XXXX", "XXX-" } ), aln_a_b );
}

BOOST_AUTO_TEST_CASE(constructs_aln_a_b_from_scaffold_lines_with_dots) {
	BOOST_CHECK_EQUAL( alignment_of_scaffold_lines( { "XXXX", "XXX." } ), aln_a_b );
}

BOOST_AUTO_TEST_CASE(constructs_aln_a_b_from_scaffold) {
	BOOST_CHECK_EQUAL( alignment_of_scaffold( "XXXX\nXXX " ), aln_a_b );
}

BOOST_AUTO_TEST_CASE(makes_correct_scaffold_lines) {
	BOOST_REQUIRE    ( ! scaffold_lines_of_alignment( aln_a_b ).empty()         );
	BOOST_CHECK_EQUAL(   scaffold_lines_of_alignment( aln_a_b ).size(),  2      );
	BOOST_CHECK_EQUAL(   scaffold_lines_of_alignment( aln_a_b ).front(), "XXXX" );
	BOOST_CHECK_EQUAL(   scaffold_lines_of_alignment( aln_a_b ).back(),  "XXX " );
}

BOOST_AUTO_TEST_CASE(scaffold_of_alignment_works_for_aln_a_b) {
	BOOST_CHECK_EQUAL( scaffold_of_alignment( aln_a_b ), "XXXX\nXXX " );
}

BOOST_AUTO_TEST_SUITE_END()
