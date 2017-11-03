/// \file
/// \brief The string_aligner test suite

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

#include "string_aligner.hpp"

#include <boost/test/unit_test.hpp>

#include <boost/numeric/conversion/cast.hpp>

#include "alignment/dyn_prog_align/detail/matrix_plotter/gnuplot_matrix_plotter.hpp"
#include "alignment/dyn_prog_align/detail/string_aligner/benchmark_dyn_prog_string_aligner.hpp"
#include "alignment/dyn_prog_align/detail/string_aligner/gen_dyn_prog_string_aligner.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/ssap_code_dyn_prog_aligner.hpp"
#include "alignment/dyn_prog_align/std_dyn_prog_aligner.hpp"
#include "alignment/gap/gap_penalty.hpp"
#include "common/algorithm/for_n.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/pair_insertion_operator.hpp"
#include "common/size_t_literal.hpp"

#include <iostream>
#include <random>

using namespace cath;
using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::align::gap;
using namespace cath::common;
using namespace std;

using boost::numeric_cast;

namespace cath {
	namespace test {

		struct string_aligner_fixture {
		protected:
			~string_aligner_fixture() noexcept = default;

			string make_random_sequence(const size_t &,
			                            const size_t &) const;

			void check_second_no_better(const str_str_score_tpl &,
			                            const str_str_score_tpl &,
			                            const gap_penalty &);

			const char_vec AA_TEST_ALPHABET = { 'A', 'C', 'D', 'E', 'F' };

			const ssap_code_dyn_prog_aligner        the_ssap_code_dyn_prog_aligner{};
			const std_dyn_prog_aligner              the_std_dyn_prog_aligner{};

			const benchmark_dyn_prog_string_aligner benchmark_string_aligner{};
			const gen_dyn_prog_string_aligner       ssap_code_string_aligner   { the_ssap_code_dyn_prog_aligner };
			const gen_dyn_prog_string_aligner       std_dyn_prog_string_aligner{ the_std_dyn_prog_aligner       };
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
string cath::test::string_aligner_fixture::make_random_sequence(const size_t &arg_min_length, ///< TODOCUMENT
                                                                const size_t &arg_max_length  ///< TODOCUMENT
                                                                ) const {
	BOOST_REQUIRE_GE(arg_max_length, arg_min_length);

	auto       rng    = default_random_engine{ random_device{}() };
	const auto length = uniform_int_distribution<size_t>{ arg_min_length, arg_max_length }( rng );

	string new_sequence;
	new_sequence.reserve(length);
	for_n(
		length,
		[&] {
			const auto random_alphabet_index = uniform_int_distribution<size_t>{ 0, AA_TEST_ALPHABET.size() - 1 }( rng );
			new_sequence.push_back( AA_TEST_ALPHABET[ random_alphabet_index ] );
		}
	);
	return new_sequence;
}

/// \brief TODOCUMENT
void cath::test::string_aligner_fixture::check_second_no_better(const str_str_score_tpl &arg_align_result_a, ///< TODOCUMENT
                                                                const str_str_score_tpl &arg_align_result_b, ///< TODOCUMENT
                                                                const gap_penalty       &arg_gap_penalty     ///< TODOCUMENT
                                                                ) {
	const score_type rescore_a = get_score_of_aligned_sequence_strings(
		get<0>( arg_align_result_a ),
		get<1>( arg_align_result_a ),
		arg_gap_penalty
	);
	const score_type rescore_b = get_score_of_aligned_sequence_strings(
		get<0>( arg_align_result_b ),
		get<1>( arg_align_result_b ),
		arg_gap_penalty
	);
	BOOST_REQUIRE_EQUAL( rescore_a, get<2>( arg_align_result_a ) );
	BOOST_CHECK_LE(      rescore_b, rescore_a );
}

BOOST_FIXTURE_TEST_SUITE(string_aligner_test_suite, cath::test::string_aligner_fixture)

/// Problems :
///  * ssap_code : if either of the sequences has length 1, this causes a crash
///  * ssap_code : doesn't achieve 1 on "AB"   / "BC"    (match "B"        with no gaps)
///  * ssap_code : doesn't achieve 1 on "FE"   / "EF"    (match "F" or "E" with no gaps)
///  * ssap_code : doesn't achieve 2 on "FACF" / "ECEFA" (match "FA"       with no gaps)

///// \brief TODOCUMENT
/////
///// "AB", "BC" with gp 0 (or indeed higher) should be able to achieve 1
//BOOST_AUTO_TEST_CASE(ssap_code_dyn_prog_finds_match_at_edge) {
//	const str_str_pair bob = align_and_format_sequence_strings(
//		"AB",
//		"BC",
//		0
//	);
//	cerr << "bob.first  : " << bob.first  << endl;
//	cerr << "bob.second : " << bob.second << endl;
//	const str_str_pair sid = align_and_format_sequence_strings(
//		"BC",
//		"AB",
//		0
//	);
//	cerr << "sid.first  : " << sid.first  << endl;
//	cerr << "sid.second : " << sid.second << endl;
//}
//
/////// \brief Check an example where, if the gap penalty is 0, gaps might get needlessly inserted
///////        to make the ends align, even with a gap penalty of 0, gaps aren't
///////
/////// \todo Ensure everything passes this
//BOOST_AUTO_TEST_CASE(handles_needless_gap_case) {
//	const str_str_pair aligned_abb_ac = dyn_prog_align( "ABB", "AC"  );
//	const str_str_pair aligned_ac_abb = dyn_prog_align( "AC",  "ABB" );
//	BOOST_CHECK_EQUAL( "ABB", aligned_abb_ac.first  );
//	BOOST_WARN_EQUAL(  "AC-", aligned_abb_ac.second );
//	BOOST_WARN_EQUAL(  "AC-", aligned_abb_ac.first  );
//	BOOST_CHECK_EQUAL( "ABB", aligned_abb_ac.second );
//}

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(basic) {
//	const size_t num_tests         = 1000;
//	const size_t min_string_length = 2;
//	const size_t max_string_length = 10;
//	for (const size_t &test_repeat_ctr : indices( num_tests ) ) {
//
//		auto         rng    = default_random_engine{ random_device{}() };
//		const auto   gap_pen = uniform_int_distribution<size_t>{ 0, max_string_length / 2 - 1 }( rng );
//		const string str1    = make_random_sequence(min_string_length, max_string_length);
//		const string str2    = make_random_sequence(min_string_length, max_string_length);
//
////        const size_t gap_pen(1);
////        const string str1("FCEEAFFED");
////        const string str2("FFCFFC"   );
////        const size_t gap_pen(1);
////        const string str1("AADCD"      );
////        const string str2("AABBBBBBCDC");
//
//		cerr << "** Gap penalty " << gap_pen << ", minimum length: " << min_string_length << ", maximum length: " << max_string_length << "**" << endl;
//		cerr << "\t" << str1.length() << "\t: " << str1 << endl;
//		cerr << "\t" << str2.length() << "\t: " << str2 << endl;
//
////		gnuplot_matrix_plotter plotter(str1.length(), str2.length(), get_window_width_for_full_matrix(str1.length(), str2.length()));
////		return_path_matrix rpm(str1.length(), str2.length(), get_window_width_for_full_matrix(str1.length(), str2.length()));
////		plotter.plot_scores( sequence_string_dyn_prog_score_source( str1, str2 ) );
////		plotter.plot_return_path_matrix( rpm );
////		plotter.finish();
//
//		const gap_penalty open_gap_pen( gap_pen, 0       );
//		const gap_penalty  all_gap_pen( gap_pen, gap_pen );
//		const str_str_score_tpl ssap_code_aligned         = ssap_code_string_aligner.align(    str1, str2, open_gap_pen );
//		const str_str_score_tpl benchmark_aligned         = benchmark_string_aligner.align(    str1, str2,  all_gap_pen );
//		const str_str_score_tpl std_dyn_prog_aligned_open = std_dyn_prog_string_aligner.align( str1, str2, open_gap_pen );
//		const str_str_score_tpl std_dyn_prog_aligned_all  = std_dyn_prog_string_aligner.align( str1, str2,  all_gap_pen );
//
////		cerr << "\tssap_code :" << endl;
////		cerr << "\t\t: " << ssap_code_aligned.get<0>()     << endl;
////		cerr << "\t\t: " << ssap_code_aligned.get<1>()     << endl;
////		cerr << "\t\t: " << ssap_code_aligned.get<2>()     << endl;
//
////		cerr << "\tbenchmark :" << endl;
////		cerr << "\t\t: " << benchmark_aligned.get<0>()     << endl;
////		cerr << "\t\t: " << benchmark_aligned.get<1>()     << endl;
////		cerr << "\t\t: " << benchmark_aligned.get<2>()     << endl;
//
////		cerr << "\tstd_dyn_prog_aligned_open :" << endl;
////		cerr << "\t\t: " << std_dyn_prog_aligned_open.get<0>() << endl;
////		cerr << "\t\t: " << std_dyn_prog_aligned_open.get<1>() << endl;
////		cerr << "\t\t: " << std_dyn_prog_aligned_open.get<2>() << endl;
//
////		cerr << "\tstd_dyn_prog_aligned_all :" << endl;
////		cerr << "\t\t: " << std_dyn_prog_aligned_all.get<0>() << endl;
////		cerr << "\t\t: " << std_dyn_prog_aligned_all.get<1>() << endl;
////		cerr << "\t\t: " << std_dyn_prog_aligned_all.get<2>() << endl;
//
////		check_second_no_better( ssap_code_aligned, benchmark_aligned, open_gap_pen );
////		check_second_no_better( benchmark_aligned, ssap_code_aligned, all_gap_pen  );
//		check_second_no_better( std_dyn_prog_aligned_open, ssap_code_aligned, open_gap_pen );
//		check_second_no_better( std_dyn_prog_aligned_all,  benchmark_aligned, all_gap_pen  );
////		throw "";
//	}
//}

/// \brief Check that get_num_gaps_and_extensions() works as expected on a few examples
BOOST_AUTO_TEST_CASE(get_num_gaps_and_extensions_works) {
	BOOST_CHECK_EQUAL( size_size_pair(1_z, 0_z), get_num_gaps_and_extensions( "A-B"       ) );
	BOOST_CHECK_EQUAL( size_size_pair(1_z, 2_z), get_num_gaps_and_extensions( "A---B"     ) );
	BOOST_CHECK_EQUAL( size_size_pair(2_z, 0_z), get_num_gaps_and_extensions( "A-B-C"     ) );
	BOOST_CHECK_EQUAL( size_size_pair(2_z, 4_z), get_num_gaps_and_extensions( "A---B---C" ) );
}

BOOST_AUTO_TEST_SUITE_END()
