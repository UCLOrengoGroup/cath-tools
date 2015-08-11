/// \file
/// \brief The alignment_breaks test suite

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

#include <boost/filesystem.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "alignment/alignment.h"
#include "alignment/io/align_scaffold.h"
#include "alignment/tools/alignment_breaks.h"
#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/type_aliases.h"
//#include "exception/invalid_argument_exception.h"
//#include "test/alignment_fixture.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;

namespace cath {
	namespace test {

		/// \brief The alignment_breaks_test_suite_fixture to assist in testing alignment_breaks
		struct alignment_breaks_test_suite_fixture {
		protected:
			~alignment_breaks_test_suite_fixture() noexcept = default;

			/// \brief An example alignment for testing alignment breaks code
			const alignment the_alignment = alignment_of_scaffold_lines( {
				// 0         1          < These two lines count up from 01 to 17 (one less than the alignment length)
				//  12345678901234567   <
				//   //     /    ////   < This line connects the between-index breaks (think: fence posts between fences) with their corresponding indices
				  "XX   XXXX  XXX   X",
				  "   XXXX  XXXX X X ",
				  "  X  XX    XX  X  "
				//   \\     \    \\\\   < This line helps re-emphasise the between-index breaks
			} );

			/// \brief The correct alignment breaks of the_alignment
			const size_vec expected_breaks{  2,  3,  9, 14, 15, 16, 17 };
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(alignment_breaks_test_suite, cath::test::alignment_breaks_test_suite_fixture)

BOOST_AUTO_TEST_CASE(alignment_breaks_are_correct_for_simple_alignment) {
	BOOST_CHECK_EQUAL_RANGES( get_alignment_breaks( the_alignment ), expected_breaks );
}

BOOST_AUTO_TEST_SUITE_END()
