/// \file
/// \brief The ssap_score_post_processing test suite

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
#include <boost/test/auto_unit_test.hpp>

#include "score/aligned_pair_score/ssap_score/ssap_score_post_processing.h"

#include "test/global_test_constants.h"

using namespace cath::score;
using namespace std;

using boost::lexical_cast;

namespace cath {
	namespace test {

		/// \brief The ssap_score_post_processing_test_suite_fixture to assist in testing ssap_score_post_processing
		struct ssap_score_post_processing_test_suite_fixture : protected global_test_constants {
		protected:
			~ssap_score_post_processing_test_suite_fixture() noexcept = default;
		};

	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(ssap_score_post_processing_test_suite, cath::test::ssap_score_post_processing_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( name_of_ssap_score_post_processing::get().size(), num_ssap_score_post_processings );
	for (const ssap_score_post_processing &value : all_ssap_score_post_processings) {
		BOOST_CHECK_GT( name_of_ssap_score_post_processing::get().at( value ).length(), 0 );
		BOOST_CHECK_GT( lexical_cast<string>( value ).length(), 0 );
	}
}

BOOST_AUTO_TEST_SUITE_END()
