/// \file
/// \brief The aligned_pair_score_list_factory test suite

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

#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/aligned_pair_score_list/aligned_pair_score_list.hpp"
#include "cath/score/aligned_pair_score_list/aligned_pair_score_list_factory.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"

using namespace ::cath::common;
using namespace ::cath::score;
using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The aligned_pair_score_list_factory_test_suite_fixture to assist in testing aligned_pair_score_list_factory
		struct aligned_pair_score_list_factory_test_suite_fixture {
		protected:
			~aligned_pair_score_list_factory_test_suite_fixture() noexcept = default;
		};

	} // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(aligned_pair_score_list_factory_test_suite, cath::test::aligned_pair_score_list_factory_test_suite_fixture)

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(short_names) {
//	const aligned_pair_score_list default_score_list = make_full_aligned_pair_score_list();
//	const str_vec                 got_short_names    = get_short_names(default_score_list);
//
//	cerr << endl;
//	for (const string &got_short_name : got_short_names) {
//		cerr << got_short_name << endl;
//	}
//	cerr << endl;
//
//	const str_vec expected_short_names = { "bob" };
//	BOOST_CHECK_EQUAL_RANGES( expected_short_names, got_short_names );
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(long_names) {
//	const aligned_pair_score_list default_score_list = make_full_aligned_pair_score_list();
//	const str_vec                 got_long_names     = get_long_names(default_score_list);
//
//	cerr << endl;
//	for (const string &got_long_name : got_long_names) {
//		cerr << got_long_name << endl;
//	}
//	cerr << endl;
//
//	const str_vec expected_long_names = { "bob" };
//	BOOST_CHECK_EQUAL_RANGES( expected_long_names, got_long_names );
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(descriptions) {
//	const aligned_pair_score_list default_score_list = make_full_aligned_pair_score_list();
//	const str_vec                 got_descriptions     = get_descriptions(default_score_list);
//
//	cerr << endl;
//	for (const string &got_description : got_descriptions) {
//		cerr << got_description << endl;
//	}
//	cerr << endl;
//
//	const str_vec expected_descriptions = { "bob" };
//	BOOST_CHECK_EQUAL_RANGES( expected_descriptions, got_descriptions );
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(references) {
//	const aligned_pair_score_list default_score_list = make_full_aligned_pair_score_list();
//	const str_vec                 got_references     = get_references(default_score_list);
//
//	cerr << endl;
//	for (const string &got_reference : got_references) {
//		cerr << got_reference << endl;
//	}
//	cerr << endl;
//
//	const str_vec expected_references = { "bob" };
//	BOOST_CHECK_EQUAL_RANGES( expected_references, got_references );
//}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
