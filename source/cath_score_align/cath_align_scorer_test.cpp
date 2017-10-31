/// \file


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

#include "cath_superpose/cath_superposer.hpp"
#include "common/argc_argv_faker.hpp"
#include "common/file/open_fstream.hpp"
#include "common/file/temp_file.hpp"
#include "common/size_t_literal.hpp"
#include "test/global_test_constants.hpp"
#include "test/predicate/files_equal.hpp"
#include "test/predicate/istream_and_file_equal.hpp"

#include <fstream>
#include <sstream>

using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The cath_align_scorer_test_suite_fixture to assist in testing cath_align_scorer
		struct cath_align_scorer_test_suite_fixture : protected global_test_constants {
		protected:
			~cath_align_scorer_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(cath_align_scorer_test_suite, cath::test::cath_align_scorer_test_suite_fixture)

/// \brief
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( 0_z, 0_z );
}

BOOST_AUTO_TEST_SUITE_END()
