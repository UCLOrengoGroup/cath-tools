/// \file
/// \brief The ssap_scores_entry test suite

#include <boost/test/auto_unit_test.hpp>

#include "common/test_tools.hpp"
#include "file/ssap_scores_file/ssap_scores_entry.hpp"

#include <vector>

using namespace cath::common;
using namespace cath::common::test;
using namespace cath::file;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The ssap_scores_entry_test_suite_fixture to assist in testing ssap_scores_entry
		struct ssap_scores_entry_test_suite_fixture {
		protected:
			~ssap_scores_entry_test_suite_fixture() noexcept = default;

			const ssap_scores_entry eg_entry{ "1cukA03", "1hjpA03", 48, 44, 94.92, 44, 91, 97, 0.71 };
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(ssap_scores_entry_test_suite, cath::test::ssap_scores_entry_test_suite_fixture)

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
