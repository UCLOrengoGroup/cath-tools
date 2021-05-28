/// \file
/// \brief The hmmer_parser test suite

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

#include <filesystem>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/resolve_hits/file/detail/hmmer_parser.hpp"
#include "cath/resolve_hits/file/parse_hmmer_out.hpp"
#include "cath/resolve_hits/full_hit_list_fns.hpp"
#include "cath/seq/seq_seg.hpp"
#include "cath/test/boost_test_print_type.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::boost::test_tools::per_element;
using ::std::filesystem::path;
using ::std::literals::string_literals::operator""s;
using ::std::make_pair;
using ::std::pair;
using ::std::string;
using ::std::vector;

namespace {

	/// \brief Type alias for a pair of string and full_hit_list
	using str_full_hit_list_pair     = pair<string, full_hit_list>;

	/// \brief Type alias for a vector of str_full_hit_list_pair values
	using str_full_hit_list_pair_vec = vector<str_full_hit_list_pair>;

	/// \brief The hmmer_parser_test_suite_fixture to assist in testing full_hit
	struct hmmer_parser_test_suite_fixture : protected global_test_constants {
	protected:
		~hmmer_parser_test_suite_fixture() noexcept = default;

		/// \brief Whether to apply CATH-specific policies
		static constexpr bool APPLY_CATH_POLICIES = false;

		/// \brief The minimum length that an alignment gap can have to be considered a gap
		static constexpr residx_t MIN_GAP_LENGTH = 3;

		/// \brief Whether to parse/output hmmsearch output alignment information
		static constexpr bool OUTPUT_HMMER_ALN = true;

		/// \brief Test whether the result of parsing the specified type of HMMER data from
		///        the specified file matches the specified data
		void test_parse(const path                       &prm_input_file,   ///< The file to parse
		                const hmmer_format               &prm_hmmer_format, ///< The type of HMMER data to expect
		                const str_full_hit_list_pair_vec &prm_expected      ///< The expected results
		                ) {
			const auto got = transform_build<str_full_hit_list_pair_vec>(
				parse_hmmer_out_file(
					prm_input_file,
					prm_hmmer_format,
					APPLY_CATH_POLICIES,
					MIN_GAP_LENGTH,
					OUTPUT_HMMER_ALN
				),
				[] (const str_calc_hit_list_pair &x) {
					return make_pair( x.first, x.second.get_full_hits() );
				}
			);

			BOOST_TEST( got == prm_expected, per_element{} );
		}

	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(hmmer_parser_test_suite, hmmer_parser_test_suite_fixture)

BOOST_AUTO_TEST_CASE(parses_jons_seqs_hmmsearch_correctly) {
	test_parse(
		CRH_HMMSCAN_DATA_DIR() / "seqs.hmmsearch",
		hmmer_format::HMMSEARCH,
		{ {
			make_pair(
				"casyuv",
				full_hit_list{ {
					full_hit{ segments_from_bounds( {   2,  93 } ), "cath|4_2_0|4xurA00/172-332-i2",  9.4, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 117, 174 } ), "cath|4_2_0|4xurA00/172-332-i2",    6, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( {   3, 120 } ), "cath|4_2_0|4xurA00/172-332-i3", 14.1, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 109, 175 } ), "cath|4_2_0|4xurA00/172-332-i3",  9.6, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( {   6, 103 } ), "cath|4_2_0|4xurA00/172-332-i4", 16.1, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 109, 175 } ), "cath|4_2_0|4xurA00/172-332-i4",  8.3, hit_score_type::BITSCORE },
				} }
			),
		} }
	);
}

BOOST_AUTO_TEST_CASE(parses_jons_seqs_hmmscan_correctly) {
	test_parse(
		CRH_HMMSCAN_DATA_DIR() / "seqs.hmmscan",
		hmmer_format::HMMSCAN,
		{ {
			make_pair(
				"casyuv",
				full_hit_list{ {
					full_hit{ segments_from_bounds( {   6, 103 } ), "cath|4_2_0|4xurA00/172-332-i4", 16.1, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 109, 175 } ), "cath|4_2_0|4xurA00/172-332-i4",  8.3, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( {   3, 120 } ), "cath|4_2_0|4xurA00/172-332-i3", 14.1, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 109, 175 } ), "cath|4_2_0|4xurA00/172-332-i3",  9.6, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( {   2,  93 } ), "cath|4_2_0|4xurA00/172-332-i2",  9.4, hit_score_type::BITSCORE },
				} }
			),
		} }
	);
}

BOOST_AUTO_TEST_CASE( detects_lack_of_alignment ) {
	try {
		parse_hmmer_out_file( CRH_HMMSCAN_DATA_DIR() / "seqs.no-alignments.hmmscan",
		                      hmmer_format::HMMSCAN,
		                      APPLY_CATH_POLICIES,
		                      MIN_GAP_LENGTH,
		                      OUTPUT_HMMER_ALN );
		BOOST_TEST( false );
	} catch ( const runtime_error_exception &ex ) {
		BOOST_TEST( ::boost::algorithm::contains( ex.what(), "Unable to find alignment data" ) );
	}
}

BOOST_AUTO_TEST_CASE(parses_jons_p53_p63_p73_hmmscan_correctly) {
	test_parse(
		CRH_HMMSCAN_DATA_DIR() / "p53_p63.hmmscan",
		hmmer_format::HMMSCAN,
		{ {
			make_pair(
				"sp|O15350|P73_HUMAN",
				full_hit_list{ {
					full_hit{ segments_from_bounds( { 111, 318 } ), "2xwcA00-i1", 327.6, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 492, 548 } ), "2eaoA01-i1",  11.7, hit_score_type::BITSCORE },
				} }
			),
			make_pair(
				"sp|P04637|P53_HUMAN",
				full_hit_list{ {
					full_hit{ segments_from_bounds( {  94, 293 } ), "3d06A00-i1", 289.4, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 324, 354 } ), "2mw4A00-i1",  16.7, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( {  28,  61 } ), "5hp0A00-i2",  12.2, hit_score_type::BITSCORE },
				} }
			),
			make_pair(
				"sp|Q9H3D4|P63_HUMAN",
				full_hit_list{ {
					full_hit{ segments_from_bounds( { 161, 367 } ), "2xwcA00-i1", 323.5, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 548, 602 } ), "2eaoA01-i2",  11.4, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 547, 601 } ), "3hilA00-i2",  11.3, hit_score_type::BITSCORE },
				} }
			)
		} }
	);
}

BOOST_AUTO_TEST_CASE(parses_ians_single_sequence_hmmscan_out_correctly) {
	test_parse(
		CRH_HMMSCAN_DATA_DIR() / "single_sequence.hmmscan.out",
		hmmer_format::HMMSCAN,
		{ {
			make_pair(
				"tr|A0A0Q0Y989|A0A0Q0Y989_9BACI",
				full_hit_list{ {
					full_hit{ segments_from_bounds( { 9,   193,           198, 217 } ), "3.20.20.10/FF/9715",  221.5, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 218, 234, 239, 333, 338, 369 } ), "2.40.37.10/FF/6607",    170, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 242,                     298 } ), "2.40.37.10/FF/6260",   25.1, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 332,                     377 } ), "2.40.37.10/FF/6260",    1.2, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 11,                      195 } ), "3.20.20.10/FF/9731",   15.9, hit_score_type::BITSCORE },
					full_hit{ segments_from_bounds( { 16,                       66 } ), "2.60.40.740/FF/2909",  10.8, hit_score_type::BITSCORE },
				} }
			),
		} }
	);
}

BOOST_AUTO_TEST_CASE(hmmsearch_input_with_hmmscan_tag) {
	BOOST_CHECK_THROW(
		test_parse(
			CRH_HMMSCAN_DATA_DIR() / "seqs.hmmsearch",
			hmmer_format::HMMSCAN,
			{}
		),
		runtime_error_exception
	);
}

BOOST_AUTO_TEST_CASE(hmmscan_input_with_hmmsearch_tag) {
	BOOST_CHECK_THROW(
		test_parse(
			CRH_HMMSCAN_DATA_DIR() / "seqs.hmmscan",
			hmmer_format::HMMSEARCH,
			{}
		),
		runtime_error_exception
	);
}

BOOST_AUTO_TEST_SUITE_END()
