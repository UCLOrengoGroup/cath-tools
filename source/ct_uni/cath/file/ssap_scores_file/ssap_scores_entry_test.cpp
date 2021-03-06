/// \file
/// \brief The ssap_scores_entry test suite

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

#include "cath/file/ssap_scores_file/ssap_scores_entry.hpp"
#include "cath/test/test_tools.hpp"

#include <vector>

using namespace ::cath::common;
using namespace ::cath::common::test;
using namespace ::cath::file;
using namespace ::std;

namespace {

	/// \brief The ssap_scores_entry_test_suite_fixture to assist in testing ssap_scores_entry
	struct ssap_scores_entry_test_suite_fixture {
	protected:
		~ssap_scores_entry_test_suite_fixture() noexcept = default;

		const ssap_scores_entry eg_entry{ "1cukA03", "1hjpA03", 48, 44, 94.92, 44, 91, 97, 0.71 };
	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(ssap_scores_entry_test_suite, ssap_scores_entry_test_suite_fixture)

BOOST_AUTO_TEST_CASE(equality_works) {
	check_equality_operators_on_diff_vals_range( vector<ssap_scores_entry>{
		eg_entry,

		ssap_scores_entry{ "white",   "1hjpA03", 48, 44, 94.92, 44, 91, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "opals",   48, 44, 94.92, 44, 91, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03",  0, 44, 94.92, 44, 91, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03", 48,  0, 94.92, 44, 91, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03", 48, 44,   0.0, 44, 91, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03", 48, 44, 94.92,  0, 91, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03", 48, 44, 94.92, 44,  0, 97, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03", 48, 44, 94.92, 44, 91,  0, 0.71 },
		ssap_scores_entry{ "1cukA03", "1hjpA03", 48, 44, 94.92, 44, 91, 97,  0.0 }
	} );
}

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL(
		to_string( eg_entry ),
		"ssap_scores_entry[1cukA03, 1hjpA03, 48, 44, 94.920000, 44, 91.000000, 97.000000, 0.710000]"
	);
}

BOOST_AUTO_TEST_CASE(parses_from_line) {
	BOOST_CHECK_EQUAL( ssap_scores_entry_from_line( "1cukA03  1hjpA03   48   44  94.92   44   91   97   0.71" ), eg_entry );
}

BOOST_AUTO_TEST_CASE(getters) {
	BOOST_CHECK_EQUAL( eg_entry.get_name_1(),     "1cukA03" );
	BOOST_CHECK_EQUAL( eg_entry.get_name_2(),     "1hjpA03" );
	BOOST_CHECK_EQUAL( eg_entry.get_length_1(),   48        );
	BOOST_CHECK_EQUAL( eg_entry.get_length_2(),   44        );
	BOOST_CHECK_EQUAL( eg_entry.get_ssap_score(), 94.92     );
	BOOST_CHECK_EQUAL( eg_entry.get_num_equivs(), 44        );
	BOOST_CHECK_EQUAL( eg_entry.get_overlap_pc(), 91        );
	BOOST_CHECK_EQUAL( eg_entry.get_seq_id_pc(),  97        );
	BOOST_CHECK_EQUAL( eg_entry.get_rmsd(),       0.71      );
}

BOOST_AUTO_TEST_SUITE_END()
