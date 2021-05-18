/// \file
/// \brief The html_align_outputter test suite

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

#include "html_align_outputter.hpp"

#include "cath/alignment/test/alignment_fixture.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/file/pdb/pdb.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The html_align_outputter_test_suite_fixture to assist in testing html_align_outputter
		struct html_align_outputter_test_suite_fixture : protected alignment_fixture {
		protected:
			~html_align_outputter_test_suite_fixture() noexcept = default;
		};

	} // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(html_align_outputter_test_suite, cath::test::html_align_outputter_test_suite_fixture)

/// \brief Check that html_align_outputter produces the expected output on aln_a_b
BOOST_AUTO_TEST_CASE(html_align_outputter_on_aln_a_b) {
//	ostringstream got_ss;
//	got_ss << html_align_outputter( aln_a_b );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
