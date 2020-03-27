/// \file
/// \brief The classn_stat test suite

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

#include "score/score_type_aliases.hpp"
#include "score/true_pos_false_neg/classn_rate_stat.hpp"
#include "score/true_pos_false_neg/true_false_pos_neg.hpp"

using namespace cath;
using namespace cath::score;

namespace cath {
	namespace test {

		/// \brief The classn_stat_test_suite_fixture to assist in testing classn_stat
		struct classn_stat_test_suite_fixture {
		protected:
			~classn_stat_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(classn_stat_test_suite, cath::test::classn_stat_test_suite_fixture)

/// \brief Check the basic classification statistics for an example true_false_pos_neg
///
/// The true_false_pos_neg has values:
///  * true  positives :  8
///  * true  negatives : 20
///  * false positives :  6
///  * false negatives : 15
///
/// Answers:
///  * sensitivity / recall / hit_rate / true_positive_rate : TP/(TP+FN) =  8 / (  8 + 15 ) =  8 / 23
///  * specificity / true_negative_rate                     : TN/(TN+FP) = 20 / ( 20 +  6 ) = 20 / 26
///  * precision / positive_predictive_value                : TP/(TP+FP) =  8 / (  8 +  6 ) =  8 / 14
///  * fall_out / false_positive_rate                       : FP/(FP+TN) =  6 / (  6 + 20 ) =  6 / 26
///  * false_discovery_rate                                 : FP/(FP+TP) =  6 / (  6 +  8 ) =  6 / 14
BOOST_AUTO_TEST_CASE(basic) {
	const true_false_pos_neg tfpn_value( 8, 20, 6, 15 );

	BOOST_CHECK_EQUAL(               sensitivity().calculate( tfpn_value ), size_rational(  8, 23 ) );
	BOOST_CHECK_EQUAL(                    recall().calculate( tfpn_value ), size_rational(  8, 23 ) );
	BOOST_CHECK_EQUAL(                  hit_rate().calculate( tfpn_value ), size_rational(  8, 23 ) );
	BOOST_CHECK_EQUAL(        true_positive_rate().calculate( tfpn_value ), size_rational(  8, 23 ) );

	BOOST_CHECK_EQUAL(               specificity().calculate( tfpn_value ), size_rational( 20, 26 ) );
	BOOST_CHECK_EQUAL(        true_negative_rate().calculate( tfpn_value ), size_rational( 20, 26 ) );

	BOOST_CHECK_EQUAL(                 precision().calculate( tfpn_value ), size_rational(  8, 14 ) );
	BOOST_CHECK_EQUAL( positive_predictive_value().calculate( tfpn_value ), size_rational(  8, 14 ) );

	BOOST_CHECK_EQUAL(                  fall_out().calculate( tfpn_value ), size_rational(  6, 26 ) );
	BOOST_CHECK_EQUAL(       false_positive_rate().calculate( tfpn_value ), size_rational(  6, 26 ) );

	BOOST_CHECK_EQUAL(      false_discovery_rate().calculate( tfpn_value ), size_rational(  6, 14 ) );
}

BOOST_AUTO_TEST_SUITE_END()

