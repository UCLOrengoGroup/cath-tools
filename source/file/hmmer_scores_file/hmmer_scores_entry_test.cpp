/// \file
/// \brief The hmmer_scores_entry test suite

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

#include "file/hmmer_scores_file/hmmer_scores_entry.hpp"
#include "test/test_tools.hpp"

#include <vector>

using namespace cath::common::test;
using namespace cath::file;
using namespace cath::file::detail;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The hmmer_scores_entry_test_suite_fixture to assist in testing hmmer_scores_entry
		struct hmmer_scores_entry_test_suite_fixture {
		protected:
			~hmmer_scores_entry_test_suite_fixture() noexcept = default;

			const hmmer_scores_entry eg_entry{ "102mA00", "-", "3ixfA00", "-", 5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-" };
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(hmmer_scores_entry_test_suite, cath::test::hmmer_scores_entry_test_suite_fixture)

BOOST_AUTO_TEST_CASE(strip_header_name_works) {
	BOOST_CHECK_EQUAL( strip_header_name( "cath|4_0_0|102mA00/0-153-i5" ), "102mA00" );
	BOOST_CHECK_EQUAL( strip_header_name( "cath|4_0_0|3ixfA00/1-137"    ), "3ixfA00" );
	BOOST_CHECK_EQUAL( strip_header_name( "cath|4_0_0|1dkgB02/138-194"  ), "1dkgB02" );
}

BOOST_AUTO_TEST_CASE(equality_works) {
	check_equality_operators_on_diff_vals_range( vector<hmmer_scores_entry>{
		eg_entry,

		hmmer_scores_entry{ "shaman",  "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "me", "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "bouquet", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "up", 5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",      0.0, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09,  0.0, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.0, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1,     0.0, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08,  0.0, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.0, 1.4, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 0.0, 1, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 0, 1, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 0, 0, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 1, 1, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 0, 1, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 0, 1, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 0, 1, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 0, "-"         },
		hmmer_scores_entry{ "102mA00", "-",  "3ixfA00", "-",  5.4e-09, 34.1, 0.1, 1.1e-08, 33.1, 0.1, 1.4, 1, 1, 0, 1, 1, 1, 1, "scattered" }
	} );
}

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_CHECK_EQUAL(
		to_string( eg_entry ),
		"hmmer_scores_entry[102mA00, -, 3ixfA00, -, 0.000000, 34.100000, 0.100000, 0.000000, 33.100000, 0.100000, 1.400000, 1, 1, 0, 1, 1, 1, 1, -]"
	);
}

BOOST_AUTO_TEST_CASE(parses_from_line) {
	BOOST_CHECK_EQUAL(
		hmmer_scores_entry_from_line( "cath|4_0_0|102mA00/0-153-i5                                   -          cath|4_0_0|3ixfA00/1-137 -            5.4e-09   34.1   0.1   1.1e-08   33.1   0.1   1.4   1   1   0   1   1   1   1 -" ),
		eg_entry
	);
}

BOOST_AUTO_TEST_CASE(getters) {
	BOOST_CHECK_EQUAL( eg_entry.get_name_1(),               "102mA00" );
	BOOST_CHECK_EQUAL( eg_entry.get_accession_1(),          "-"       );
	BOOST_CHECK_EQUAL( eg_entry.get_name_2(),               "3ixfA00" );
	BOOST_CHECK_EQUAL( eg_entry.get_accession_2(),          "-"       );
	BOOST_CHECK_EQUAL( eg_entry.get_full_sequence_evalue(), 5.4e-09   );
	BOOST_CHECK_EQUAL( eg_entry.get_full_sequence_score(),  34.1      );
	BOOST_CHECK_EQUAL( eg_entry.get_full_sequence_bias(),   0.1       );
	BOOST_CHECK_EQUAL( eg_entry.get_best_1_domain_evalue(), 1.1e-08   );
	BOOST_CHECK_EQUAL( eg_entry.get_best_1_domain_score(),  33.1      );
	BOOST_CHECK_EQUAL( eg_entry.get_best_1_domain_bias(),   0.1       );
	BOOST_CHECK_EQUAL( eg_entry.get_expected_num_doms(),    1.4       );
	BOOST_CHECK_EQUAL( eg_entry.get_reg(),                  1         );
	BOOST_CHECK_EQUAL( eg_entry.get_clu(),                  1         );
	BOOST_CHECK_EQUAL( eg_entry.get_ov(),                   0         );
	BOOST_CHECK_EQUAL( eg_entry.get_env(),                  1         );
	BOOST_CHECK_EQUAL( eg_entry.get_dom(),                  1         );
	BOOST_CHECK_EQUAL( eg_entry.get_rep(),                  1         );
	BOOST_CHECK_EQUAL( eg_entry.get_inc(),                  1         );
	BOOST_CHECK_EQUAL( eg_entry.get_description(),          "-"       );
}

BOOST_AUTO_TEST_SUITE_END()

