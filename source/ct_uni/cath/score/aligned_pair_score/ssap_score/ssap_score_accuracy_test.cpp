/// \file
/// \brief The ssap_score_accuracy test suite

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
#include <boost/test/unit_test.hpp>

#include "cath/score/aligned_pair_score/ssap_score/ssap_score_accuracy.hpp"

#include "cath/test/global_test_constants.hpp"

using namespace ::cath::score;
using namespace ::std;

using ::boost::lexical_cast;

/// \brief TODOCUMENT
BOOST_AUTO_TEST_SUITE(ssap_score_accuracy_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( name_of_ssap_score_accuracy::get().size(), num_ssap_score_accuracies );
	for (const ssap_score_accuracy &value : all_ssap_score_accuracies) {
		BOOST_CHECK_GT( name_of_ssap_score_accuracy::get().at( value ).length(), 0 );
		BOOST_CHECK_GT( lexical_cast<string>( value ).length(), 0 );
	}
}

BOOST_AUTO_TEST_SUITE_END()
