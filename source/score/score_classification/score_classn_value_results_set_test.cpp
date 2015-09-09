/// \file
/// \brief The score_classn_value_results_set test suite

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

#include <boost/algorithm/string/join.hpp> // ***** TEMPORARY *****
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>
#include <boost/range/algorithm_ext/insert.hpp>
#include <boost/range/algorithm_ext/insert.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/algorithm/contains.h"
#include "common/algorithm/sort_uniq_copy.h"
#include "common/algorithm/transform_build.h"
#include "common/file/simple_file_read_write.h"
#include "file/ssap_scores_file/ssap_scores_entry.h"
#include "file/ssap_scores_file/ssap_scores_entry_to_score_classn_value.h"
#include "file/ssap_scores_file/ssap_scores_file.h"
#include "score/pair_scatter_plotter/pair_scatter_plotter.h"
#include "score/score_classification/label_pair_is_positive/label_pair_is_positive.h"
#include "score/true_pos_false_neg/classn_rate_stat.h"
#include "score/true_pos_false_neg/classn_stat_pair_series.h"
#include "score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter.h"
#include "score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.h"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list.h"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list_list.h"
#include "score_classn_value_results_set.h"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

using boost::algorithm::icontains;
using boost::algorithm::join; // ***** TEMPORARY *****
using boost::filesystem::path;
using boost::none;
using boost::optional;
using boost::range::insert;
using cath::file::ssap_scores_entry;
using cath::file::ssap_scores_file;

namespace cath {
	namespace test {

		/// \brief The score_classn_value_results_set_test_suite_fixture to assist in testing score_classn_value_results_set
		struct score_classn_value_results_set_test_suite_fixture {
		protected:
			~score_classn_value_results_set_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			static bool is_positive(const string &arg_string
			                        ) {
				static str_vec positives;
				if ( positives.empty() ) {
					positives = sort_copy( read_file<string>( "/cath/homes2/ucbctnl/ssap_invest_over_subset_data/random_pairs_subset.positives.txt" ) );
				}
				return contains( positives, arg_string );
			}

			/// \brief TODOCUMENT
			const score_classn_value_results_set & get_som2_scores() const {
				static score_classn_value_results_set some2_scores;
				if ( some2_scores.empty() ) {
					some2_scores = read_from_dir( path( "/export/people/ucbctnl/som2_scores" ), select_file() );
				}
				return some2_scores;
			}

			/// \brief TODOCUMENT
			const score_classn_value_results_set & get_refined_som2_scores() const {
				static score_classn_value_results_set refined_som2_scores;
				if ( refined_som2_scores.empty() ) {
					refined_som2_scores = read_from_dir( path( "/export/people/ucbctnl/refined_som2_scores" ), select_file() );
				}
				return refined_som2_scores;
			}

			void plot_both(const score_classn_value_list_vec    &arg_results,        ///< TODOCUMENT
			               const string                         &arg_file_stem_stem, ///< TODOCUMENT
			               const vector<pair<string, opt_str>>  &arg_series_to_plot  ///< TODOCUMENT
			               ) {
				plot_roc(
					classn_stat_plotter(),
					path( arg_file_stem_stem + ".roc" ),
					arg_results,
					make_standard_score_roc_plotter_spec( arg_series_to_plot )
				);

				plot_precision_recall(
					classn_stat_plotter(),
					path( arg_file_stem_stem + ".precision_recall" ),
					arg_results,
					make_standard_score_precision_recall_plotter_spec( arg_series_to_plot )
				);
			}

			/// \brief TODOCUMENT
			class select_file {
			public:
				optional<pair<bool, string>> operator()(const path &arg_file) const {
//					optional<pair<bool, string>> select_file(const path &arg_file) {
					using opt_bool_str_pair = optional<pair<bool, string>>;
					if ( file_size( arg_file ) < 100
							|| icontains( arg_file.string(), "stderr" )
							|| icontains( arg_file.string(), "svmlight_data" ) ) {
						return opt_bool_str_pair( none );
					}

					const auto file_stem = arg_file.stem().string();
					return opt_bool_str_pair{ make_pair(
						score_classn_value_results_set_test_suite_fixture::is_positive( file_stem ),
						file_stem
					) };
				}
			};
		};

		/// \brief TODOCUMENT
		struct ticket_913_fixture : protected score_classn_value_results_set_test_suite_fixture {
			const path root_dir                 = path( "/cath/homes2/ucbctnl/svm_gubbins_files" );
//			const path random_pairs_subset_file = root_dir / "random_pairs_subset.txt";

			const path labelled_pair_list       = root_dir / "pair_list.labelled";

			const path hmmscan_results_file     = root_dir / "results.hmmscan";
			const path hmmsearch_results_file   = root_dir / "results.hmmsearch";
			const path comp_prc_results_dir     = root_dir / "results.prc";
			const path ssap_results_file        = root_dir / "results.ssap";



			// grep -PR '\d+\s+cath\|' . -H | awk '{print $1 " " $4 " " $2}' | tr '|' ' ' | tr ':' ' ' | tr '/' ' ' | sed 's/\.out//g' | awk '{print $2 " " $7 " " $4}' | grep -Fwf $ESU/svm/svm_data/random_pairs_subset.txt > $ESU/svm/svm_data/hhsearch_comp_results/hhsearch_comp_full_results.txt
			//
			// ls prc_results_with_duplicates_removed/* | xargs awk '{print $1 " " $6 " BREAK " $0}' | grep -Fwf random_pairs_subset.txt | sort | cut -d ' ' -f 4- > prc_results_with_duplicates_removed/prc_results_with_duplicates_removed.txt
			// awk '{print $1 " " $2 " BREAK " $0}' ssap_results/ssap_results | grep -Fwf random_pairs_subset.txt | sort | cut -d ' ' -f 4- > ssap_results.sorted.txt
			//
			// cd /cath/godzilla-data1/people/ucbctnl/hhsuite_pairwise_v4_0_0_benchmark/johannes_results/
			// grep -PR '\d+\s+cath\|' . -H | sed 's/^\.\///g' | sed 's/\.out:/ /g' | perl -n -e 'my $line = $_; chomp($line); if ( $line =~ /^\b(\d\w{3}\S\d{2})\b.*\b(\d\w{3}\S\d{2})\b/ ) { print "$1 $2 BREAK $line\n"; } else { die $line }' | grep -Fwf /export/people/ucbctnl/svm/svm_data/random_pairs_subset.sorted.txt | sort | cut -d ' ' -f 4- > /export/people/ucbctnl/svm/svm_data/jo
			//
			// sort hmmscan_results/hmmscan_results.txt     > hmmscan_results/hmmscan_results.sorted.txt
			// sort hmmsearch_results/hmmsearch_results.txt > hmmsearch_results/hmmsearch_results.sorted.txt

			// Prefer hmmscan over hmmsearch

			// Build magic function from SSAP and PRC
			// SSAP score alone
			// SSAP score plus gubbins
		};

	}
}



/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(score_classn_value_results_set_test_suite, cath::test::score_classn_value_results_set_test_suite_fixture)




BOOST_FIXTURE_TEST_SUITE(ticket_913, cath::test::ticket_913_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ssap) {
	BOOST_CHECK( true );
//	const auto ssap_data = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file );
//	BOOST_REQUIRE_EQUAL( ssap_data.size(), 21576 );

//	const auto is_pos = make_label_pair_is_positive( labelled_pair_list );
//
//	const score_classn_value_list_vec the_results = {
//		make_val_list_of_ssap_scores_entries( ssap_data, is_pos, &ssap_scores_entry::get_ssap_score, true,  "SSAP score" ),
//		make_val_list_of_ssap_scores_entries( ssap_data, is_pos, &ssap_scores_entry::get_rmsd,       false, "RMSD"       )
//	};
//
//	plot_both(
//		the_results,
//		"/tmp/saxophone_right_here",
//		{
//			{ "SSAP score", opt_str{ R"( linetype 1 linecolor rgb "black" linewidth 3 )" } },
//			{ "RMSD",       opt_str{ R"( linetype 1 linecolor rgb "red"   linewidth 3 )" } }
//		}
//	);
}

BOOST_AUTO_TEST_SUITE_END()





//BOOST_AUTO_TEST_SUITE(score_classn_value_results_set_test_suite)
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/1.ssap",
//		{ { "ssap", opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } } }
//	);
//
////	plot_roc(
////		classn_stat_plotter(),
////		"/export/people/ucbctnl/sab_graphs/1.ssap.roc",
////		som2_results,
////		make_standard_score_roc_plotter_spec( { { "ssap", opt_str{ R"( linetype 1 linecolor rgb "#000000" )" } } } )
////	);
////
////	plot_precision_recall(
////		classn_stat_plotter(),
////		"/export/people/ucbctnl/sab_graphs/1.ssap.precision_recall",
////		som2_results,
////		make_standard_score_precision_recall_plotter_spec( { { "ssap", opt_str{ R"( linetype 1 linecolor rgb "#000000" )" } } } )
////	);
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_naives) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/2.ssap_plus_naives",
//		{
//			{ "ssap",                                                 opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "shorter_protein_length",                               opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "geometric_mean_protein_length",                        opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length", opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                  opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "longer_protein_length",                                opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } }
//		}
//	);
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_basics) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/3.ssap_plus_basics",
//		{
//			{ "ssap",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "GSAS",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "SAS",                                                                                    opt_str{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "SIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               opt_str{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "num_aligned_residues",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     opt_str{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    opt_str{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//		}
//	);
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_natives_and_non_natives) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//
//	const score_classn_value_list native_tm_align_results = read_score_classn_value_list(
//		"/export/people/ucbctnl/tm_align_data.txt",
//		true,
//		"Native TM-align",
//		[&] (const str_vec &x) { return is_positive( x.front() ); }
//	);
//
//	som2_results.push_back( native_tm_align_results );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/4.ssap_plus_natives_and_non_natives",
//		{
//			{ "ssap",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
////			{ "GSAS",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
////			{ "SAS",                                                                                    opt_str{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#2222ff" linewidth 3 )" } },
////			{ "MIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
////			{ "SIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
////			{ "sequence_id.identity",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#006401"             )" } },
////			{ "lDDT.threshold[STD_MEAN]",                                                               opt_str{ R"( linetype 1 linecolor rgb "#010067"             )" } },
////			{ "num_aligned_residues",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
////			{ "overlap.num_aligned_residues.longer_protein_length",                                     opt_str{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
////			{ "dRMSD",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
////			{ "shorter_protein_length",                                                                 opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
////			{ "SI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
////			{ "MI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
////			{ "RMSD",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
////			{ "geometric_mean_protein_length",                                                          opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
////			{ "overlap.shorter_protein_length.longer_protein_length",                                   opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
////			{ "mean_protein_length",                                                                    opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
////			{ "sequence_id.blosum62",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
////			{ "longer_protein_length",                                                                  opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
////			{ "overlap.num_aligned_residues.shorter_protein_length",                                    opt_str{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//			{ "Native TM-align",                                                                        opt_str{ R"( linetype 1 linecolor rgb "#000088" linewidth 3 )" } }
//		}
//	);
//}
//
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_natives_and_non_natives_on_refined) {
//	auto orig_results = get_refined_som2_scores();
//	auto som2_results = make_score_classn_value_list_vec( orig_results );
//
//	const score_classn_value_list native_tm_align_results = read_score_classn_value_list(
//		"/export/people/ucbctnl/tm_align_data.txt",
//		true,
//		"Native TM-align",
//		[&] (const str_vec &x) { return is_positive( x.front() ); }
//	);
//
//	som2_results.push_back( native_tm_align_results );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/5.ssap_plus_natives_and_non_natives_on_refined",
//		{
//			{ "ssap",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
////			{ "GSAS",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
////			{ "SAS",                                                                                    opt_str{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#2222ff" linewidth 3 )" } },
////			{ "MIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
////			{ "SIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
////			{ "sequence_id.identity",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#006401"             )" } },
////			{ "lDDT.threshold[STD_MEAN]",                                                               opt_str{ R"( linetype 1 linecolor rgb "#010067"             )" } },
////			{ "num_aligned_residues",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
////			{ "overlap.num_aligned_residues.longer_protein_length",                                     opt_str{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
////			{ "dRMSD",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
////			{ "shorter_protein_length",                                                                 opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
////			{ "SI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
////			{ "MI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
////			{ "RMSD",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
////			{ "geometric_mean_protein_length",                                                          opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
////			{ "overlap.shorter_protein_length.longer_protein_length",                                   opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
////			{ "mean_protein_length",                                                                    opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
////			{ "sequence_id.blosum62",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
////			{ "longer_protein_length",                                                                  opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
////			{ "overlap.num_aligned_residues.shorter_protein_length",                                    opt_str{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//			{ "Native TM-align",                                                                        opt_str{ R"( linetype 1 linecolor rgb "#000088" linewidth 3 )" } }
//		}
//	);
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_basics_plus_extras) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/6.ssap_plus_basics_plus_extras",
//		{
//			{ "ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified", opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70].cb_atoms",                                            opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "ssap",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "SAS.select_best_score_percent[70].cb_atoms",                                             opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70].cb_atoms",                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70]",                                                 opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70]",                                                     opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_best_score_percent[70]",                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.cb_atoms",                                                                          opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01].cb_atoms",                                                   opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01].cb_atoms",                                               opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01]",                                                            opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "structal.cb_atoms",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.select_min_score[0.01].cb_atoms",                                               opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.cb_atoms",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01].cb_atoms",                                                    opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01]",                                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.cb_atoms",                                                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01]",                                                             opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS",                                                                                    opt_str{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score.select_min_score[0.01]",                                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX.select_best_score_percent[70].cb_atoms",                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX.select_best_score_percent[70].cb_atoms",                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "MIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "lDDT.threshold[4]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[2]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               opt_str{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "lDDT.threshold[1]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "num_aligned_residues",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "lDDT.threshold[0.5]",                                                                    opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     opt_str{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    opt_str{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//		}
//	);
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_basics_plus_extras_plus_svm) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//
//	const auto svm_light_play = read_svmlight_predictions_files( { {
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.1.test_predictions"  ),  "SVM01"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.2.test_predictions"  ),  "SVM02"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.3.test_predictions"  ),  "SVM03"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.4.test_predictions"  ),  "SVM04"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.5.test_predictions"  ),  "SVM05"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.6.test_predictions"  ),  "SVM06"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.7.test_predictions"  ),  "SVM07"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.8.test_predictions"  ),  "SVM08"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.9.test_predictions"  ),  "SVM09"  },
//		{  path(  "/export/people/ucbctnl/som2_scores/some_svmlight_data.10.test_predictions" ),  "SVM10"  }
//	} } );
//
//	insert( som2_results, end( som2_results ), svm_light_play );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/7.ssap_plus_basics_plus_extras_plus_svm",
//		{
//			{ "SVM01",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM02",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM03",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM04",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM05",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM06",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM07",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM08",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM09",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM10",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//
//
//			{ "ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified", opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70].cb_atoms",                                            opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "ssap",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "SAS.select_best_score_percent[70].cb_atoms",                                             opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70].cb_atoms",                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70]",                                                 opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70]",                                                     opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_best_score_percent[70]",                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.cb_atoms",                                                                          opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01].cb_atoms",                                                   opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01].cb_atoms",                                               opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01]",                                                            opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "structal.cb_atoms",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.select_min_score[0.01].cb_atoms",                                               opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.cb_atoms",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01].cb_atoms",                                                    opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01]",                                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.cb_atoms",                                                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01]",                                                             opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS",                                                                                    opt_str{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score.select_min_score[0.01]",                                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX.select_best_score_percent[70].cb_atoms",                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX.select_best_score_percent[70].cb_atoms",                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "MIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "lDDT.threshold[4]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[2]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               opt_str{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "lDDT.threshold[1]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "num_aligned_residues",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "lDDT.threshold[0.5]",                                                                    opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     opt_str{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    opt_str{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//		}
//	);
//}
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap_plus_basics_plus_extras_plus_svm_on_refined) {
//	auto orig_results = get_refined_som2_scores();
//	auto som2_results = make_score_classn_value_list_vec( orig_results );
//
//	const auto svm_light_play = read_svmlight_predictions_files( { {
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.1.test_predictions"  ),  "SVM01"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.2.test_predictions"  ),  "SVM02"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.3.test_predictions"  ),  "SVM03"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.4.test_predictions"  ),  "SVM04"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.5.test_predictions"  ),  "SVM05"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.6.test_predictions"  ),  "SVM06"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.7.test_predictions"  ),  "SVM07"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.8.test_predictions"  ),  "SVM08"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.9.test_predictions"  ),  "SVM09"  },
//		{  path(  "/export/people/ucbctnl/refined_som2_scores/svmlight_data.10.test_predictions" ),  "SVM10"  }
//	} } );
////
//	insert( som2_results, end( som2_results ), svm_light_play );
//
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/8.ssap_plus_basics_plus_extras_plus_svm_on_refined",
//		{
//			{ "SVM01",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM02",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM03",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM04",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM05",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM06",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM07",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM08",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM09",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM10",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//
//
//			{ "ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified", opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70].cb_atoms",                                            opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "ssap",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "SAS.select_best_score_percent[70].cb_atoms",                                             opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70].cb_atoms",                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70]",                                                 opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70]",                                                     opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_best_score_percent[70]",                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.cb_atoms",                                                                          opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01].cb_atoms",                                                   opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01].cb_atoms",                                               opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01]",                                                            opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "structal.cb_atoms",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.select_min_score[0.01].cb_atoms",                                               opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.cb_atoms",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01].cb_atoms",                                                    opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01]",                                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.cb_atoms",                                                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01]",                                                             opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS",                                                                                    opt_str{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score.select_min_score[0.01]",                                                        opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score",                                                                               opt_str{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX.select_best_score_percent[70].cb_atoms",                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX.select_best_score_percent[70].cb_atoms",                                           opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "MIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "lDDT.threshold[4]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[2]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               opt_str{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "lDDT.threshold[1]",                                                                      opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "num_aligned_residues",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "lDDT.threshold[0.5]",                                                                    opt_str{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     opt_str{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  opt_str{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 opt_str{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     opt_str{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   opt_str{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          opt_str{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   opt_str{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    opt_str{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   opt_str{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  opt_str{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    opt_str{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//		}
//	);
//
////	mt19937 rng{ random_device{}() };
////	write_to_svm_light_data_files( orig_results, "/export/people/ucbctnl/refined_som2_scores/svmlight_data", 10, rng, 0.5 );
//}


///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(dssp_based_ssap_roc_and_pr) {
//	// cat /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_alignments/*.stdout | awk '{print $1 "_" $2 " " $5}' > /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt
//	// /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt
//
//	const score_classn_value_list dssp_ssap_results = read_score_classn_value_list(
//		"/cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt",
//		true,
//		"dssp_based_cath-ssap",
//		[&] (const str_vec &x) { return is_positive( x.front() ); }
//	);
//
//	const score_classn_value_list r17937_ssap_results = read_score_classn_value_list(
//		"/cath/homes2/ucbctnl/ssap_invest_over_subset_data/old_ssap_simple_summary.score_classn_value_list.txt",
//		true,
//		"SSAP_r17937",
//		[&] (const str_vec &x) { return is_positive( x.front() ); }
//	);
//
//	const score_classn_value_list r7906_ssap_results = read_score_classn_value_list(
//		"/cath/homes2/ucbctnl/ssap_invest_over_subset_data/r7906_ssap_simple_summary.score_classn_value_list.txt",
//		true,
//		"SSAP_r7906",
//		[&] (const str_vec &x) { return is_positive( x.front() ); }
//	);
//
//	auto som2_results = get_som2_scores();
//
//	const bool warn_on_missing = true;
//	add_score_classn_value_list_and_add_missing(
//		som2_results,
//		dssp_ssap_results,
//		worst_possible_score( dssp_ssap_results ),
//		warn_on_missing
//	);
//	add_score_classn_value_list_and_add_missing(
//		som2_results,
//		r17937_ssap_results,
//		worst_possible_score( r17937_ssap_results ),
//		warn_on_missing
//	);
//	add_score_classn_value_list_and_add_missing(
//		som2_results,
//		r7906_ssap_results,
//		worst_possible_score( r7906_ssap_results ),
//		warn_on_missing
//	);
//
//	plot_both(
//		make_score_classn_value_list_vec( som2_results ),
//		"/tmp/dssp_discrimination",
//		{
//			{ "ssap",                 opt_str{ R"( linetype 1 linecolor rgb "red"    linewidth 3 )" } },
//			{ "dssp_based_cath-ssap", opt_str{ R"( linetype 1 linecolor rgb "blue"   linewidth 3 )" } },
//			{ "SSAP_r17937",          opt_str{ R"( linetype 1 linecolor rgb "green"  linewidth 3 )" } },
//			{ "SSAP_r7906",           opt_str{ R"( linetype 1 linecolor rgb "orange" linewidth 3 )" } }
//		}
//	);
//}

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(dssp_based_ssap_score_correlation) {
//	// cat /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_alignments/*.stdout | awk '{print $1 "_" $2 " " $5}' > /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt
//	// /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt
//
//	const score_classn_value_list dssp_ssap_results = read_score_classn_value_list(
//		"/cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt",
//		true,
//		"dssp_based_cath-ssap",
//		[&] (const str_vec &x) { return is_positive( x.front() ); }
//	);
//
//	auto som2_results = get_som2_scores();
//
//	const bool warn_on_missing = true;
//	add_score_classn_value_list_and_add_missing(
//		som2_results,
//		dssp_ssap_results,
//		worst_possible_score( dssp_ssap_results ),
//		warn_on_missing
//	);
//
//	const auto the_correl_data = get_correlated_data( som2_results, "ssap", "dssp_based_cath-ssap" );
//	pair_scatter_plotter().plot( "/tmp/dssp_scatter_plot", the_correl_data, "ssap", "dssp_based_cath-ssap" );
//}


///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(basic) {
//	/// FN should take a single path and return optional<pair<bool, string>>,
//	/// representing files that should be processed with a bool indicating whether it's a positive instance and a label
//
//	// sed -i '/Missing some residues whilst loading alignment/d' *.stderr
//	// sed -i '/Scores are as follows/d' *.stderr
//
//	// cat /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_alignments/*.stdout | awk '{print $1 "_" $2 " " $5}' | grep -vP '1a0jA01_1gpzA02|1arbA02_1gpzA02|1aym400_1vq8L01|1befA02_1gpzA02|1eb7A01_1tme400|1fod400_2mev400|1gpzA02_1k28D01|1txkA01_1vq8L01' > /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt
//	// /cath/homes2/ucbctnl/ssap_invest_over_subset_data/dssp_simple_summary.score_classn_value_list.txt
//
////	get_instance_labels();
//
//
////	auto lovely_results = read_from_dir(
//////		path( "/cath/homes2/ucbctnl/ssap_invest_over_subset_data/many_scores" ),
//////		path( "/export/people/ucbctnl/many_scores" ),
//////		path( "/export/people/ucbctnl/full_scores" ),
////		path( "/export/people/ucbctnl/som2_scores" ),
////		[&] (const path &arg_file) {
////			using opt_bool_str_pair = optional<pair<bool, string>>;
////			if ( file_size( arg_file ) < 100 || icontains( arg_file.string(), "stderr" ) ) {
////				return opt_bool_str_pair( none );
////			}
////
////			const auto file_stem = arg_file.stem().string();
////			return opt_bool_str_pair{ make_pair(
////				is_positive( file_stem ),
////				file_stem
////			) };
////		}
////	);
////	mt19937 rng{ random_device{}() };
////	write_to_svm_light_data_files( lovely_results, "/tmp/some_svmlight_data", 10, rng, 0.5 );
//
//
//
//
//
//////	const auto &ssap__match = find_score_classn_value_list_of_name( lovely_results, "ssap" );
//////	summarise_score_classn_value_list( cerr, ssap__match );
////
////	const auto named_true_false_pos_neg_lists = make_named_true_false_pos_neg_list_list( lovely_results );
////
////	if ( named_true_false_pos_neg_lists.empty() ) {
////		BOOST_THROW_EXCEPTION(runtime_error_exception("Failed to read any scores files for score_classn_value_results_set_test"));
////	}
////
////	const auto aucs = sort_copy(
////		areas_under_roc_curves( named_true_false_pos_neg_lists ),
////		[] (const str_doub_pair &x, const str_doub_pair &y) {
////			return ( x.second < y.second );
////		}
////	);
////
////	for (const auto &auc : aucs) {
////		cerr << "AUC : "<< std::right << std::setw( 10 ) << auc.second << " " << auc.first << endl;
////	}
//////	cerr << join( aucs | lexical_casted<string>(), ", " ) << endl;
//
//
//	const auto svm_light_play = read_svmlight_predictions_files( { {
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.0001_c0.0001220703125.test_predictions"     ),  "RBF  [gamma:  0.0001;  c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.0001_c0.015625.test_predictions"            ),  "RBF  [gamma:  0.0001;  c:  c0.015625"             },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.0001_c0.25.test_predictions"                ),  "RBF  [gamma:  0.0001;  c:  c0.25"                 },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.0001_c32.test_predictions"                  ),  "RBF  [gamma:  0.0001;  c:  c32"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.001_c0.000244140625.test_predictions"       ),  "RBF  [gamma:  0.001;   c:  c0.000244140625"       },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.001_c0.00048828125.test_predictions"        ),  "RBF  [gamma:  0.001;   c:  c0.00048828125"        },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.001_c0.015625.test_predictions"             ),  "RBF  [gamma:  0.001;   c:  c0.015625"             },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.001_c7.62939453125e-06.test_predictions"    ),  "RBF  [gamma:  0.001;   c:  c7.62939453125e-06"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.01_c0.000244140625.test_predictions"        ),  "RBF  [gamma:  0.01;    c:  c0.000244140625"       },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.01_c0.0009765625.test_predictions"          ),  "RBF  [gamma:  0.01;    c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.01_c0.00390625.test_predictions"            ),  "RBF  [gamma:  0.01;    c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.01_c0.125.test_predictions"                 ),  "RBF  [gamma:  0.01;    c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.01_c7.62939453125e-06.test_predictions"     ),  "RBF  [gamma:  0.01;    c:  c7.62939453125e-06"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.1_c0.0001220703125.test_predictions"        ),  "RBF  [gamma:  0.1;     c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.1_c0.000244140625.test_predictions"         ),  "RBF  [gamma:  0.1;     c:  c0.000244140625"       },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.1_c0.0078125.test_predictions"              ),  "RBF  [gamma:  0.1;     c:  c0.0078125"            },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.1_c0.03125.test_predictions"                ),  "RBF  [gamma:  0.1;     c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_0.1_c0.5.test_predictions"                    ),  "RBF  [gamma:  0.1;     c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100000_c0.0001220703125.test_predictions"     ),  "RBF  [gamma:  100000;  c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100000_c0.0078125.test_predictions"           ),  "RBF  [gamma:  100000;  c:  c0.0078125"            },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100000_c0.0625.test_predictions"              ),  "RBF  [gamma:  100000;  c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100000_c16.test_predictions"                  ),  "RBF  [gamma:  100000;  c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100000_c32.test_predictions"                  ),  "RBF  [gamma:  100000;  c:  c32"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100000_c8.test_predictions"                   ),  "RBF  [gamma:  100000;  c:  c8"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10000_c0.0001220703125.test_predictions"      ),  "RBF  [gamma:  10000;   c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10000_c0.0009765625.test_predictions"         ),  "RBF  [gamma:  10000;   c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10000_c1.52587890625e-05.test_predictions"    ),  "RBF  [gamma:  10000;   c:  c1.52587890625e-05"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10000_c3.0517578125e-05.test_predictions"     ),  "RBF  [gamma:  10000;   c:  c3.0517578125e-05"     },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10000_c6.103515625e-05.test_predictions"      ),  "RBF  [gamma:  10000;   c:  c6.103515625e-05"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10000_c8.test_predictions"                    ),  "RBF  [gamma:  10000;   c:  c8"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1000_c0.0009765625.test_predictions"          ),  "RBF  [gamma:  1000;    c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1000_c0.03125.test_predictions"               ),  "RBF  [gamma:  1000;    c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1000_c3.814697265625e-06.test_predictions"    ),  "RBF  [gamma:  1000;    c:  c3.814697265625e-06"   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1000_c6.103515625e-05.test_predictions"       ),  "RBF  [gamma:  1000;    c:  c6.103515625e-05"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1000_c7.62939453125e-06.test_predictions"     ),  "RBF  [gamma:  1000;    c:  c7.62939453125e-06"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100_c0.001953125.test_predictions"            ),  "RBF  [gamma:  100;     c:  c0.001953125"          },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100_c0.00390625.test_predictions"             ),  "RBF  [gamma:  100;     c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_100_c0.0625.test_predictions"                 ),  "RBF  [gamma:  100;     c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c0.0001220703125.test_predictions"         ),  "RBF  [gamma:  10;      c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c0.0009765625.test_predictions"            ),  "RBF  [gamma:  10;      c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c0.00390625.test_predictions"              ),  "RBF  [gamma:  10;      c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c0.0625.test_predictions"                  ),  "RBF  [gamma:  10;      c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c0.125.test_predictions"                   ),  "RBF  [gamma:  10;      c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c1.9073486328125e-06.test_predictions"     ),  "RBF  [gamma:  10;      c:  c1.9073486328125e-06"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c6.103515625e-05.test_predictions"         ),  "RBF  [gamma:  10;      c:  c6.103515625e-05"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_10_c7.62939453125e-06.test_predictions"       ),  "RBF  [gamma:  10;      c:  c7.62939453125e-06"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1_c0.0078125.test_predictions"                ),  "RBF  [gamma:  1;       c:  c0.0078125"            },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1_c0.03125.test_predictions"                  ),  "RBF  [gamma:  1;       c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1_c0.5.test_predictions"                      ),  "RBF  [gamma:  1;       c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c0.000244140625.test_predictions"       ),  "RBF  [gamma:  1e-05;   c:  c0.000244140625"       },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c0.03125.test_predictions"              ),  "RBF  [gamma:  1e-05;   c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c0.125.test_predictions"                ),  "RBF  [gamma:  1e-05;   c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c0.25.test_predictions"                 ),  "RBF  [gamma:  1e-05;   c:  c0.25"                 },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c0.5.test_predictions"                  ),  "RBF  [gamma:  1e-05;   c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c1.test_predictions"                    ),  "RBF  [gamma:  1e-05;   c:  c1"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-05_c8.test_predictions"                    ),  "RBF  [gamma:  1e-05;   c:  c8"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-06_c1.9073486328125e-06.test_predictions"  ),  "RBF  [gamma:  1e-06;   c:  c1.9073486328125e-06"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-07_c0.125.test_predictions"                ),  "RBF  [gamma:  1e-07;   c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-07_c4.test_predictions"                    ),  "RBF  [gamma:  1e-07;   c:  c4"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c0.0001220703125.test_predictions"      ),  "RBF  [gamma:  1e-08;   c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c0.00390625.test_predictions"           ),  "RBF  [gamma:  1e-08;   c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c0.03125.test_predictions"              ),  "RBF  [gamma:  1e-08;   c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c0.0625.test_predictions"               ),  "RBF  [gamma:  1e-08;   c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c0.125.test_predictions"                ),  "RBF  [gamma:  1e-08;   c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c1.52587890625e-05.test_predictions"    ),  "RBF  [gamma:  1e-08;   c:  c1.52587890625e-05"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c1.9073486328125e-06.test_predictions"  ),  "RBF  [gamma:  1e-08;   c:  c1.9073486328125e-06"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-08_c4.test_predictions"                    ),  "RBF  [gamma:  1e-08;   c:  c4"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-09_c0.00390625.test_predictions"           ),  "RBF  [gamma:  1e-09;   c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-09_c0.0625.test_predictions"               ),  "RBF  [gamma:  1e-09;   c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-10_c0.000244140625.test_predictions"       ),  "RBF  [gamma:  1e-10;   c:  c0.000244140625"       },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-10_c0.25.test_predictions"                 ),  "RBF  [gamma:  1e-10;   c:  c0.25"                 },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-10_c0.5.test_predictions"                  ),  "RBF  [gamma:  1e-10;   c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-10_c16.test_predictions"                   ),  "RBF  [gamma:  1e-10;   c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-10_c4.test_predictions"                    ),  "RBF  [gamma:  1e-10;   c:  c4"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-10_c9.5367431640625e-07.test_predictions"  ),  "RBF  [gamma:  1e-10;   c:  c9.5367431640625e-07"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c0.0001220703125.test_predictions"      ),  "RBF  [gamma:  1e-11;   c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c0.001953125.test_predictions"          ),  "RBF  [gamma:  1e-11;   c:  c0.001953125"          },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c0.00390625.test_predictions"           ),  "RBF  [gamma:  1e-11;   c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c0.015625.test_predictions"             ),  "RBF  [gamma:  1e-11;   c:  c0.015625"             },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c0.125.test_predictions"                ),  "RBF  [gamma:  1e-11;   c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c1.9073486328125e-06.test_predictions"  ),  "RBF  [gamma:  1e-11;   c:  c1.9073486328125e-06"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-11_c3.0517578125e-05.test_predictions"     ),  "RBF  [gamma:  1e-11;   c:  c3.0517578125e-05"     },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c0.0001220703125.test_predictions"      ),  "RBF  [gamma:  1e-12;   c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c0.0078125.test_predictions"            ),  "RBF  [gamma:  1e-12;   c:  c0.0078125"            },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c0.015625.test_predictions"             ),  "RBF  [gamma:  1e-12;   c:  c0.015625"             },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c0.125.test_predictions"                ),  "RBF  [gamma:  1e-12;   c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c0.5.test_predictions"                  ),  "RBF  [gamma:  1e-12;   c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c16.test_predictions"                   ),  "RBF  [gamma:  1e-12;   c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-12_c3.0517578125e-05.test_predictions"     ),  "RBF  [gamma:  1e-12;   c:  c3.0517578125e-05"     },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-13_c0.00390625.test_predictions"           ),  "RBF  [gamma:  1e-13;   c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-13_c0.0625.test_predictions"               ),  "RBF  [gamma:  1e-13;   c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-13_c0.5.test_predictions"                  ),  "RBF  [gamma:  1e-13;   c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-13_c1.9073486328125e-06.test_predictions"  ),  "RBF  [gamma:  1e-13;   c:  c1.9073486328125e-06"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-13_c6.103515625e-05.test_predictions"      ),  "RBF  [gamma:  1e-13;   c:  c6.103515625e-05"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c0.5.test_predictions"                  ),  "RBF  [gamma:  1e-14;   c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c16.test_predictions"                   ),  "RBF  [gamma:  1e-14;   c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c32.test_predictions"                   ),  "RBF  [gamma:  1e-14;   c:  c32"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c4.test_predictions"                    ),  "RBF  [gamma:  1e-14;   c:  c4"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c6.103515625e-05.test_predictions"      ),  "RBF  [gamma:  1e-14;   c:  c6.103515625e-05"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c8.test_predictions"                    ),  "RBF  [gamma:  1e-14;   c:  c8"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-14_c9.5367431640625e-07.test_predictions"  ),  "RBF  [gamma:  1e-14;   c:  c9.5367431640625e-07"  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c0.000244140625.test_predictions"       ),  "RBF  [gamma:  1e-15;   c:  c0.000244140625"       },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c0.0009765625.test_predictions"         ),  "RBF  [gamma:  1e-15;   c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c0.001953125.test_predictions"          ),  "RBF  [gamma:  1e-15;   c:  c0.001953125"          },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c0.015625.test_predictions"             ),  "RBF  [gamma:  1e-15;   c:  c0.015625"             },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c0.0625.test_predictions"               ),  "RBF  [gamma:  1e-15;   c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c1.test_predictions"                    ),  "RBF  [gamma:  1e-15;   c:  c1"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c2.test_predictions"                    ),  "RBF  [gamma:  1e-15;   c:  c2"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-15_c32.test_predictions"                   ),  "RBF  [gamma:  1e-15;   c:  c32"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-16_c0.001953125.test_predictions"          ),  "RBF  [gamma:  1e-16;   c:  c0.001953125"          },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-16_c0.03125.test_predictions"              ),  "RBF  [gamma:  1e-16;   c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-16_c0.5.test_predictions"                  ),  "RBF  [gamma:  1e-16;   c:  c0.5"                  },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-17_c0.0078125.test_predictions"            ),  "RBF  [gamma:  1e-17;   c:  c0.0078125"            },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-17_c0.03125.test_predictions"              ),  "RBF  [gamma:  1e-17;   c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-17_c3.814697265625e-06.test_predictions"   ),  "RBF  [gamma:  1e-17;   c:  c3.814697265625e-06"   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-17_c6.103515625e-05.test_predictions"      ),  "RBF  [gamma:  1e-17;   c:  c6.103515625e-05"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c0.0001220703125.test_predictions"      ),  "RBF  [gamma:  1e-18;   c:  c0.0001220703125"      },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c0.03125.test_predictions"              ),  "RBF  [gamma:  1e-18;   c:  c0.03125"              },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c0.0625.test_predictions"               ),  "RBF  [gamma:  1e-18;   c:  c0.0625"               },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c0.25.test_predictions"                 ),  "RBF  [gamma:  1e-18;   c:  c0.25"                 },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c16.test_predictions"                   ),  "RBF  [gamma:  1e-18;   c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c7.62939453125e-06.test_predictions"    ),  "RBF  [gamma:  1e-18;   c:  c7.62939453125e-06"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-18_c8.test_predictions"                    ),  "RBF  [gamma:  1e-18;   c:  c8"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-19_c0.00048828125.test_predictions"        ),  "RBF  [gamma:  1e-19;   c:  c0.00048828125"        },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-19_c0.0009765625.test_predictions"         ),  "RBF  [gamma:  1e-19;   c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-19_c0.125.test_predictions"                ),  "RBF  [gamma:  1e-19;   c:  c0.125"                },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-19_c16.test_predictions"                   ),  "RBF  [gamma:  1e-19;   c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-19_c1.test_predictions"                    ),  "RBF  [gamma:  1e-19;   c:  c1"                    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-19_c7.62939453125e-06.test_predictions"    ),  "RBF  [gamma:  1e-19;   c:  c7.62939453125e-06"    },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-20_c0.0009765625.test_predictions"         ),  "RBF  [gamma:  1e-20;   c:  c0.0009765625"         },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-20_c0.00390625.test_predictions"           ),  "RBF  [gamma:  1e-20;   c:  c0.00390625"           },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-20_c16.test_predictions"                   ),  "RBF  [gamma:  1e-20;   c:  c16"                   },
////		{  path(  "/export/people/ucbctnl/some_svmlight/some_svmlight_data.2.rbf_gamma_1e-20_c1.9073486328125e-06.test_predictions"  ),  "RBF  [gamma:  1e-20;   c:  c1.9073486328125e-06"  }
//
//		{  path(  "/tmp/some_svmlight_data.1.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.2.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.3.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.4.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.5.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.6.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.7.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.8.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.9.test_predictions"  ),  "SVM [RBF;gamma:1;c:5]"  },
//		{  path(  "/tmp/some_svmlight_data.10.test_predictions" ),  "SVM [RBF;gamma:1;c:5]"  }
//	} } );
//
//	const auto svm_light_aucs = sort_copy(
//		areas_under_roc_curves( make_named_true_false_pos_neg_list_list( svm_light_play ) ),
//		[] (const str_doub_pair &x, const str_doub_pair &y) {
//			return ( x.second < y.second );
//		}
//	);
//
//	for (const auto &svm_light_auc : svm_light_aucs) {
//		cerr << "AUC : "<< std::right << std::setw( 10 ) << svm_light_auc.second << " " << svm_light_auc.first << endl;
//	}
//
////	const auto svmlight_play                 = read_svmlight_predictions_file( "/tmp/lovely_svm_data.1.test_predictions", "svmlight_play" );
////	const auto svmlight_play_auc             = area_under_roc_curve( svmlight_play );
////	const auto svmlight_play_named_tfpn_list = make_named_true_false_pos_neg_list( svmlight_play );
////	cerr << "Area under svmlight_play is " << svmlight_play_auc << endl;
////	const auto roc_series_list = classn_stat_pair_series( svmlight_play_named_tfpn_list );
//	plot_roc(
//		classn_stat_plotter(),
//		"/tmp/lovely_svm_plot",
//		svm_light_play,
//		make_standard_score_roc_plotter_spec()
//	);
//
////	plot_precision_recall(
////		classn_stat_plotter(),
////		"/tmp/lovely_prec_rec_plot",
////		svm_light_play,
////		{}
////	);
//
//
//
////	plot_roc(
////		classn_stat_plotter(),
////		"/tmp/some_svmlight_plot",
////		lovely_results,
////		{}
////	);
//}
//
//BOOST_AUTO_TEST_SUITE_END()





//AUC :   0.546841 overlap.num_aligned_residues.shorter_protein_length
//AUC :   0.605014 longer_protein_length
//AUC :   0.613516 sequence_id.blosum62
//AUC :   0.655347 mean_protein_length
//AUC :   0.662545 overlap.shorter_protein_length.longer_protein_length
//AUC :   0.670037 geometric_mean_protein_length
//AUC :   0.702592 RMSD
//AUC :   0.704105 MI
//AUC :   0.706146 SI
//AUC :   0.709255 shorter_protein_length
//AUC :   0.709281 dRMSD
//AUC :   0.712885 overlap.num_aligned_residues.longer_protein_length
//AUC :   0.713382 lDDT.threshold[0.5]
//AUC :   0.726469 num_aligned_residues
//AUC :   0.746452 lDDT.threshold[1]
//AUC :   0.758863 lDDT.threshold[STD_MEAN]
//AUC :   0.764754 lDDT.threshold[2]
//AUC :   0.765439 sequence_id.identity
//AUC :   0.767857 SIMAX
//AUC :   0.769322 lDDT.threshold[4]
//AUC :   0.771355 MIMAX
//AUC :   0.784733 SIMAX.select_best_score_percent[70].cb_atoms
//AUC :   0.789956 MIMAX.select_best_score_percent[70].cb_atoms
//AUC :   0.794369 TM-score
//AUC :   0.794589 TM-score.select_min_score[0.01]
//AUC :   0.796118 structal
//AUC :   0.796122 SAS
//AUC :   0.796318 SAS.select_min_score[0.01]
//AUC :   0.796733 SAS.cb_atoms
//AUC :   0.796785 structal.select_min_score[0.01]
//AUC :   0.796805 SAS.select_min_score[0.01].cb_atoms
//AUC :   0.796821 TM-score.cb_atoms
//AUC :   0.797073 TM-score.select_min_score[0.01].cb_atoms
//AUC :    0.79756 structal.cb_atoms
//AUC :   0.798055 GSAS
//AUC :   0.798063 GSAS.select_min_score[0.01]
//AUC :   0.798211 structal.select_min_score[0.01].cb_atoms
//AUC :   0.798774 GSAS.select_min_score[0.01].cb_atoms
//AUC :   0.798958 GSAS.cb_atoms
//AUC :   0.809053 SAS.select_best_score_percent[70]
//AUC :   0.810459 GSAS.select_best_score_percent[70]
//AUC :   0.810802 structal.select_best_score_percent[70]
//AUC :   0.812547 structal.select_best_score_percent[70].cb_atoms
//AUC :   0.813306 SAS.select_best_score_percent[70].cb_atoms
//AUC :   0.814712 ssap
//AUC :    0.81486 GSAS.select_best_score_percent[70].cb_atoms
//AUC :   0.822689 ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified

BOOST_AUTO_TEST_SUITE_END()

