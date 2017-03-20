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
#include <boost/log/trivial.hpp>
#include <boost/optional.hpp>
#include <boost/optional.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm_ext/insert.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/irange.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/set_union_build.hpp"
#include "common/algorithm/sort_uniq_copy.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/file/simple_file_read_write.hpp"
#include "file/hmmer_scores_file/hmmer_scores_entry.hpp"
#include "file/hmmer_scores_file/hmmer_scores_file.hpp"
#include "file/prc_scores_file/prc_scores_entry.hpp"
#include "file/prc_scores_file/prc_scores_file.hpp"
#include "file/ssap_scores_file/ssap_scores_entry.hpp"
#include "file/ssap_scores_file/ssap_scores_entry_to_score_classn_value.hpp"
#include "file/ssap_scores_file/ssap_scores_file.hpp"
#include "score/pair_scatter_plotter/pair_scatter_plotter.hpp"
#include "score/score_classification/label_pair_is_positive/label_pair_is_positive.hpp"
#include "score/score_classification/value_list_scaling.hpp"
#include "score/true_pos_false_neg/classn_rate_stat.hpp"
#include "score/true_pos_false_neg/classn_stat_pair_series.hpp"
#include "score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter.hpp"
#include "score/true_pos_false_neg/classn_stat_plotter/classn_stat_plotter_spec.hpp"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list.hpp"
#include "score/true_pos_false_neg/named_true_false_pos_neg_list_list.hpp"
#include "score_classn_value_results_set.hpp"

#include <cmath>
#include <functional>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::score;
using namespace std;

using boost::adaptors::map_keys;
using boost::algorithm::icontains;
using boost::algorithm::join; // ***** TEMPORARY *****
using boost::filesystem::path;
using boost::irange;
using boost::is_space;
using boost::make_optional;
using boost::none;
using boost::numeric_cast;
using boost::optional;
using boost::range::combine;
using boost::range::copy;
using boost::range::insert;
using boost::token_compress_on;

namespace cath {
	namespace test {

		/// \brief The score_classn_value_results_set_test_suite_fixture to assist in testing score_classn_value_results_set
		struct score_classn_value_results_set_test_suite_fixture {
		protected:
			~score_classn_value_results_set_test_suite_fixture() noexcept = default;

			/// \brief TODOCUMENT
			static bool is_positive(const string &arg_string ///< TODOCUMENT
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
			               const vector<pair<string, str_opt>>  &arg_series_to_plot  ///< TODOCUMENT
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
					using bool_str_pair_opt = optional<pair<bool, string>>;
					if ( file_size( arg_file ) < 100
							|| icontains( arg_file.string(), "stderr" )
							|| icontains( arg_file.string(), "svmlight_data" ) ) {
						return bool_str_pair_opt( none );
					}

					const auto file_stem = arg_file.stem().string();
					return bool_str_pair_opt{ make_pair(
						score_classn_value_results_set_test_suite_fixture::is_positive( file_stem ),
						file_stem
					) };
				}
			};
		};

		/// \brief TODOCUMENT
		struct ticket_913_fixture : protected score_classn_value_results_set_test_suite_fixture {
			const path root_dir                 = path( "/cath/homes2/ucbctnl/svm_gubbins_files" );
//			const path root_dir                 = path( "/home/lewis/svm_gubbins_files" );
			const path random_pairs_subset_file = root_dir / "random_pairs_subset.txt";

			const path svm_data_dir             = root_dir / "svm_experiments_data";

			const path graphs_dir               = root_dir / "graphs";

			const path labelled_pair_list       = root_dir / "pair_list.labelled";

			const path hmmscan_results_file     = root_dir / "results.hmmscan";
			const path hmmsearch_results_file   = root_dir / "results.hmmsearch";
			const path comp_prc_results_file    = root_dir / "results.prc";
			const path ssap_results_file        = root_dir / "results.ssap";
			const path multi_struc_results_file = root_dir / "results.multi_struc_scores";

			const path ssap_svm_data_dir                = svm_data_dir / "ssap";
			const path prc_svm_data_dir                 = svm_data_dir / "prc";
			const path hmmer_svm_data_dir               = svm_data_dir / "hmmer";
			const path ssap_and_prc_svm_data_dir        = svm_data_dir / "ssap_and_prc";
			const path ssap_and_hmmer_svm_data_dir        = svm_data_dir / "ssap_and_hmmer";
			const path ssap_prc_and_hmmer_svm_data_dir  = svm_data_dir / "ssap_prc_and_hmmer";
			const path ssap_prc_and_multi_svm_data_dir  = svm_data_dir / "ssap_prc_and_multi";

			const path ssap_svm_data_stem               = ssap_svm_data_dir               / "svm_data";
			const path prc_svm_data_stem                = prc_svm_data_dir                / "svm_data";
			const path hmmer_svm_data_stem              = hmmer_svm_data_dir              / "svm_data";
			const path ssap_and_prc_svm_data_stem       = ssap_and_prc_svm_data_dir       / "svm_data";
			const path ssap_and_hmmer_svm_data_stem     = ssap_and_hmmer_svm_data_dir     / "svm_data";
			const path ssap_prc_and_hmmer_svm_data_stem = ssap_prc_and_hmmer_svm_data_dir / "svm_data";
			const path ssap_prc_and_multi_svm_data_stem = ssap_prc_and_multi_svm_data_dir / "svm_data";

			static auto normal_label_getter() {
				return [] (const auto &x) {
					return make_pair( x.get_name_1(), x.get_name_2() );
				};
			}

			/// \brief TODOCUMENT
			static str_vec_vec parse_multi_struc_scores(const path &arg_filename ///< TODOCUMENT
			                                            ) {
				ifstream input_stream;
				common::open_ifstream( input_stream, arg_filename );
				string line_string;
				str_vec_vec results;
				while ( getline( input_stream, line_string ) ) {
					results.push_back(
						split_build<str_vec>( line_string, is_space(), token_compress_on )
					);
				}
				input_stream.close();
				return results;
			}

			/// \brief TODOCUMENT
			///
			/// \todo Ideally, make this a variadic template
			static score_classn_value_list_vec join_value_lists(const score_classn_value_list_vec &arg_list_a, ///< TODOCUMENT
			                                                    const score_classn_value_list_vec &arg_list_b  ///< TODOCUMENT
			                                                    ) {
				score_classn_value_list_vec results;
				results.reserve( arg_list_a.size() + arg_list_b.size() );
				copy( arg_list_a, back_inserter( results ) );
				copy( arg_list_b, back_inserter( results ) );
				return results;
			}

			/// \brief TODOCUMENT
			template <typename... Ts>
			static score_classn_value_list_vec join_value_lists(const score_classn_value_list_vec &arg_list_a, ///< TODOCUMENT
			                                                    const score_classn_value_list_vec &arg_list_b, ///< TODOCUMENT
			                                                    const Ts &...                      arg_lists   ///< TODOCUMENT
			                                                    ) {
				return join_value_lists( join_value_lists( arg_list_a, arg_list_b ), arg_lists... );
			}

			/// \brief TODOCUMENT
			static score_classn_value_list_vec get_multi_ssap_data(const ssap_scores_entry_vec  &arg_ssap_data, ///< TODOCUMENT
			                                                       const label_pair_is_positive &arg_is_pos     ///< TODOCUMENT
			                                                       ) {
//				const auto name_getter = [] (const ssap_scores_entry &x) { return make_pair( x.get_name_1(), x.get_name_2() ); };
				const auto name_getter = normal_label_getter();
				return {
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, true,  "ssap.length1"    , [] (const ssap_scores_entry  &x) { return x.get_length_1  ();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, true,  "ssap.length2"    , [] (const ssap_scores_entry  &x) { return x.get_length_2  ();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, true,  "ssap.ssap_score" , [] (const ssap_scores_entry  &x) { return x.get_ssap_score();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, true,  "ssap.num_equivs" , [] (const ssap_scores_entry  &x) { return x.get_num_equivs();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, true,  "ssap.overlap_pc" , [] (const ssap_scores_entry  &x) { return x.get_overlap_pc();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, true,  "ssap.seq_id_pc"  , [] (const ssap_scores_entry  &x) { return x.get_seq_id_pc ();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data, arg_is_pos, false, "ssap.rmsd"       , [] (const ssap_scores_entry  &x) { return x.get_rmsd      ();      }, name_getter )
				};
			}

			/// \brief TODOCUMENT
			static score_classn_value_list_vec value_lists_of_ssap_data(const ssap_scores_entry_vec  &arg_ssap_data, ///< TODOCUMENT
			                                                            const label_pair_is_positive &arg_is_pos     ///< TODOCUMENT
			                                                            ) {
//				const auto name_getter = [] (const ssap_scores_entry &x) { return make_pair( x.get_name_1(), x.get_name_2() ); };
				const auto name_getter = normal_label_getter();
				return {
//					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, true,  "ssap.length1"    , [] (const ssap_scores_entry  &x) { return x.get_length_1  ();      }, name_getter ),
//					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, true,  "ssap.length2"    , [] (const ssap_scores_entry  &x) { return x.get_length_2  ();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, true,  "ssap.ssap_score" , [] (const ssap_scores_entry  &x) { return x.get_ssap_score();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, true,  "ssap.num_equivs" , [] (const ssap_scores_entry  &x) { return x.get_num_equivs();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, true,  "ssap.overlap_pc" , [] (const ssap_scores_entry  &x) { return x.get_overlap_pc();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, true,  "ssap.seq_id_pc"  , [] (const ssap_scores_entry  &x) { return x.get_seq_id_pc ();      }, name_getter ),
					make_val_list_of_scores_entries( arg_ssap_data,  arg_is_pos, false, "ssap.rmsd"       , [] (const ssap_scores_entry  &x) { return x.get_rmsd      ();      }, name_getter )
				};
			}

			/// \brief TODOCUMENT
			static score_classn_value_list_vec value_lists_of_prc_data(const prc_scores_entry_vec   &arg_prc_data, ///< TODOCUMENT
			                                                           const label_pair_is_positive &arg_is_pos    ///< TODOCUMENT
			                                                           ) {
//				const auto name_getter = [] (const prc_scores_entry &x) { return make_pair( x.get_name_1(), x.get_name_2() ); };
				const auto name_getter = normal_label_getter();
				return {
//					make_val_list_of_scores_entries( arg_prc_data,   arg_is_pos, true,  "prc.length1"     , [] (const prc_scores_entry   &x) { return        x.get_length_1(); }, name_getter ),
//					make_val_list_of_scores_entries( arg_prc_data,   arg_is_pos, true,  "prc.length2"     , [] (const prc_scores_entry   &x) { return        x.get_length_2(); }, name_getter ),
					make_val_list_of_scores_entries( arg_prc_data,   arg_is_pos, false, "prc.evalue"      , [] (const prc_scores_entry   &x) { return log10( x.get_evalue() ); }, name_getter ),
					make_val_list_of_scores_entries( arg_prc_data,   arg_is_pos, true,  "prc.simple"      , [] (const prc_scores_entry   &x) { return        x.get_simple();   }, name_getter ),
					make_val_list_of_scores_entries( arg_prc_data,   arg_is_pos, true,  "prc.reverse"     , [] (const prc_scores_entry   &x) { return        x.get_reverse();  }, name_getter )

				};
			}

			/// \brief TODOCUMENT
			static score_classn_value_list_vec value_lists_of_hmmer_data(const hmmer_scores_entry_vec &arg_hmmer_data, ///< TODOCUMENT
			                                                             const label_pair_is_positive &arg_is_pos      ///< TODOCUMENT
			                                                             ) {
				const auto name_getter = normal_label_getter();
				return {
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, false, "hmmer.full_sequence_evalue", [] (const hmmer_scores_entry &x) { return log10( x.get_full_sequence_evalue() ); }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.full_sequence_score",  [] (const hmmer_scores_entry &x) { return        x.get_full_sequence_score();    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.full_sequence_bias",   [] (const hmmer_scores_entry &x) { return        x.get_full_sequence_bias();     }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, false, "hmmer.best_1_domain_evalue", [] (const hmmer_scores_entry &x) { return log10( x.get_best_1_domain_evalue() ); }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.best_1_domain_score",  [] (const hmmer_scores_entry &x) { return        x.get_best_1_domain_score();    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.best_1_domain_bias",   [] (const hmmer_scores_entry &x) { return        x.get_best_1_domain_bias();     }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.expected_num_doms",    [] (const hmmer_scores_entry &x) { return        x.get_expected_num_doms();      }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.reg",                  [] (const hmmer_scores_entry &x) { return        x.get_reg();                    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.clu",                  [] (const hmmer_scores_entry &x) { return        x.get_clu();                    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.ov",                   [] (const hmmer_scores_entry &x) { return        x.get_ov();                     }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.env",                  [] (const hmmer_scores_entry &x) { return        x.get_env();                    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.dom",                  [] (const hmmer_scores_entry &x) { return        x.get_dom();                    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.rep",                  [] (const hmmer_scores_entry &x) { return        x.get_rep();                    }, name_getter ),
					make_val_list_of_scores_entries( arg_hmmer_data, arg_is_pos, true,  "hmmer.inc",                  [] (const hmmer_scores_entry &x) { return        x.get_inc();                    }, name_getter )
				};
			}

			/// \brief TODOCUMENT
			static score_classn_value_list_vec value_lists_of_multi_data(const str_vec_vec            &arg_multi_data, ///< TODOCUMENT
																	     const label_pair_is_positive &arg_is_pos      ///< TODOCUMENT
																	     ) {
				return transform_build<score_classn_value_list_vec>(
					irange( 2_z, arg_multi_data.front().size() ), // Start from 2 to skip the two IDs
					[&] (const size_t &idx) {
						return make_val_list_of_scores_entries(
							arg_multi_data,
							arg_is_pos,
							true,
							"multi_struc_score_" + to_string( idx - 1 ),
							[&] (const str_vec &row) { return stod( row[ idx ] ); },
							[&] (const str_vec &row) { return make_pair( row[ 0 ], row[ 1 ] ); }
						);
					}
				);
			}

			/// \brief TODOCUMENT
			static score_classn_value_list_vec magic_function_data(const ssap_scores_entry_vec  &arg_ssap_data, ///< TODOCUMENT,
			                                                       const prc_scores_entry_vec   &arg_prc_data,  ///< TODOCUMENT
			                                                       const label_pair_is_positive &arg_is_pos     ///< TODOCUMENT
			                                                       ) {
//				map<str_str_pair, reference_wrapper<const ssap_scores_entry>> ssaps_by_label;
//				for (const auto &x : arg_ssap_data) {
//					ssaps_by_label.emplace( make_pair( x.get_name_1(), x.get_name_2() ), std::cref( x ) );
//				}
//				map<str_str_pair, reference_wrapper<const prc_scores_entry>> prcs_by_label;
//				for (const auto &x : arg_prc_data) {
//					prcs_by_label.emplace( make_pair( x.get_name_1(), x.get_name_2() ), std::cref( x ) );
//				}
				const auto ssaps_by_label = transform_build<map<str_str_pair, reference_wrapper<const ssap_scores_entry>>>(
					arg_ssap_data,
					[] (const ssap_scores_entry &x) {
						return make_pair( make_pair( x.get_name_1(), x.get_name_2() ), std::cref( x ) );
					}
				);
				const auto prcs_by_label = transform_build<map<str_str_pair, reference_wrapper<const prc_scores_entry>>>(
					arg_prc_data,
					[] (const prc_scores_entry &x) {
						return make_pair( make_pair( x.get_name_1(), x.get_name_2() ), std::cref( x ) );
					}
				);
				const auto labels = set_union_build<set<str_str_pair>>(
					ssaps_by_label | map_keys,
					prcs_by_label  | map_keys
				);
				const auto magic_function_results = make_val_list_of_scores_entries(
					labels,
					arg_is_pos,
					true,
					"magic.function",
					[&] (const str_str_pair &x) {
						const double &ssap_score = ssaps_by_label.at( x ).get().get_ssap_score();
						const auto    prc_itr    = prcs_by_label.find( x );
						if ( prc_itr == common::cend( prcs_by_label ) ) {
							return numeric_limits<double>::lowest();
						};
						const double &prc_score = prc_itr->second.get().get_evalue();
						return ssap_score - log10( prc_score );
					},
					[&] (const str_str_pair &x) {
						const auto &ssap_result = ssaps_by_label.at( x ).get();
						return make_pair( ssap_result.get_name_1(), ssap_result.get_name_2() );
					}
				);
				return { { magic_function_results } };
			}

			/// \brief TODOCUMENT
			///
			/// \todo Ideally, make this a variadic template
			static vector<pair<string, str_opt>> join_series_lists(const vector<pair<string, str_opt>> &arg_list_a, ///< TODOCUMENT
			                                                       const vector<pair<string, str_opt>> &arg_list_b  ///< TODOCUMENT
			                                                       ) {
				vector<pair<string, str_opt>> results;
				results.reserve( arg_list_a.size() + arg_list_b.size() );
				copy( arg_list_a, back_inserter( results ) );
				copy( arg_list_b, back_inserter( results ) );
				return results;
			}

			/// \brief TODOCUMENT
			template <typename... Ts>
			static vector<pair<string, str_opt>> join_series_lists(const vector<pair<string, str_opt>> &arg_list_a, ///< TODOCUMENT
			                                                       const vector<pair<string, str_opt>> &arg_list_b, ///< TODOCUMENT
			                                                       const Ts &...                        arg_lists   ///< TODOCUMENT
			                                                       ) {
				return join_series_lists( join_series_lists( arg_list_a, arg_list_b ), arg_lists... );
			}

//			/// \brief TODOCUMENT
//			///
//			/// \todo Ideally, make this a variadic template
//			static vector<pair<string, str_opt>> join_series_lists(const vector<pair<string, str_opt>> &arg_list_a, ///< TODOCUMENT
//			                                                       const vector<pair<string, str_opt>> &arg_list_b, ///< TODOCUMENT
//			                                                       const vector<pair<string, str_opt>> &arg_list_c  ///< TODOCUMENT
//			                                                       ) {
//				return join_series_lists( join_series_lists( arg_list_a, arg_list_b ), arg_list_c );
//			}
//
//			/// \brief TODOCUMENT
//			///
//			/// \todo Ideally, make this a variadic template
//			static vector<pair<string, str_opt>> join_series_lists(const vector<pair<string, str_opt>> &arg_list_a, ///< TODOCUMENT
//			                                                       const vector<pair<string, str_opt>> &arg_list_b, ///< TODOCUMENT
//			                                                       const vector<pair<string, str_opt>> &arg_list_c, ///< TODOCUMENT
//			                                                       const vector<pair<string, str_opt>> &arg_list_d  ///< TODOCUMENT
//			                                                       ) {
//				return join_series_lists( join_series_lists( join_series_lists( arg_list_a, arg_list_b ), arg_list_c ), arg_list_d );
//			}

			// "black"   linewidth 3
			// "#00ff00"
			// "#0000ff"
			// "#ff0000"
			// "#01fffe"
			// "#ffa6fe"
			// "#ffdb66"
			// "#006401"
			// "#010067"
			// "#95003a"
			// "#007db5"
			// "#ff00f6"
			// "#774d00"
			// "#90fb92"
			// "#0076ff"
			// "#d5ff00"
			// "#ff937e"
			// "#6a826c"
			// "#ff029d"
			// "#fe8900"
			// "#7a4782"
			// "#7e2dd2"

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> ssap_series_list() {
				return {
//					{ "ssap.length1",    str_opt{ R"( linetype 1 linecolor rgb "black"               )" } },
//					{ "ssap.length2",    str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
					{ "ssap.ssap_score", str_opt{ R"( linetype 1 linecolor rgb "#0000ff" linewidth 3 )" } },
					{ "ssap.num_equivs", str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
					{ "ssap.overlap_pc", str_opt{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
					{ "ssap.seq_id_pc",  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
					{ "ssap.rmsd",       str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } }
				};
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> prc_series_list() {
				return {
//					{ "prc.length1", str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//					{ "prc.length2", str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
					{ "prc.evalue",  str_opt{ R"( linetype 1 linecolor rgb "#95003a" linewidth 3 )" } },
					{ "prc.simple",  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "prc.reverse", str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } }
				};
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> hmmer_series_list() {
				return {
					{ "hmmer.full_sequence_evalue", str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.full_sequence_score",  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.full_sequence_bias",   str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.best_1_domain_evalue", str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.best_1_domain_score",  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.best_1_domain_bias",   str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.expected_num_doms",    str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.reg",                  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.clu",                  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.ov",                   str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.env",                  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.dom",                  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.rep",                  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
					{ "hmmer.inc",                  str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } }
				};
			}


			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> magic_series_list() {
				return {
					{ "magic.function", str_opt{ R"( linetype 1 linecolor rgb "#774d00" linewidth 3 )" } }
				};
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> make_svm_series_list(const string &arg_name_stem,   ///< TODOCUMENT
			                                                          const string &arg_spec_string, ///< TODOCUMENT
			                                                          const size_t &arg_num          ///< TODOCUMENT
			                                                          ) {
				return transform_build<vector<pair<string, str_opt>>>(
					irange( 0_z, arg_num ),
					[&] (const size_t &x) {
						const string title_string = x == 0 ? string{ ""          }
						                                   : string{ " notitle " };
						return make_pair(
							arg_name_stem + to_string( x + 1 ),
							make_optional( R"( linetype 1 linewidth 1 )" + arg_spec_string + title_string )
						);
					}
				);
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__prc_alone_series_list() {
				return make_svm_series_list( "svm__prc_alone", R"( linecolor rgb "#FFCC00" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__ssap_alone_series_list() {
				return make_svm_series_list( "svm__ssap_alone", R"( linecolor rgb "#FF029D" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__hmmer_alone_series_list() {
				return make_svm_series_list( "svm__hmmer_alone", R"( linecolor rgb "#FE8900" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__ssap_prc_series_list() {
				return make_svm_series_list( "svm__ssap_prc", R"( linecolor rgb "#22FF22" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__ssap_hmmer_series_list() {
				return make_svm_series_list( "svm__ssap_hmmer", R"( linecolor rgb "#01fffe" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__ssap_prc_hmmer_series_list() {
				return make_svm_series_list( "svm__ssap_prc_hmmer", R"( linecolor rgb "#2222FF" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__ssap_prc_multi_series_list() {
				return make_svm_series_list( "svm__ssap_prc_multi", R"( linecolor rgb "#FF7F00" )", 5 );
			}

			/// \brief TODOCUMENT
			static vector<pair<string, str_opt>> svm__full_series_list() {
				return {
					{ "svm__full", str_opt{ R"( linetype 1 linewidth 4 linecolor rgb "#000000" )" } }
				};
			}

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

	}  // namespace test
}  // namespace cath



/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(score_classn_value_results_set_test_suite, cath::test::score_classn_value_results_set_test_suite_fixture)

BOOST_FIXTURE_TEST_SUITE(ticket_913, cath::test::ticket_913_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK( true );
}

//BOOST_AUTO_TEST_CASE(make_data_from_ssap) {
//	const auto ssap_data    = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file );
//	const auto is_pos       = make_label_pair_is_positive( labelled_pair_list );
//	const auto full_results = value_lists_of_ssap_data   ( ssap_data, is_pos  );
//	const auto the_set      = make_score_classn_value_results_set( full_results );
//
//	auto rng = mt19937{ random_device{}() };
//	write_to_svm_light_data_files(
//		the_set,
//		ssap_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//BOOST_AUTO_TEST_CASE(make_data_from_prc) {
//	const auto prc_data     = prc_scores_file::parse_prc_scores_file_fancy( comp_prc_results_file );
//	const auto is_pos       = make_label_pair_is_positive( labelled_pair_list );
//	const auto full_results = value_lists_of_prc_data( prc_data, is_pos );
//	const auto the_set      = make_score_classn_value_results_set( full_results );
//
//	BOOST_LOG_TRIVIAL( warning ) << "The PRC results don't cover all pairs and the code doesn't yet handle that";
//
//	auto rng = mt19937{ random_device{}() };
//	write_to_svm_light_data_files(
//		the_set,
//		prc_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//BOOST_AUTO_TEST_CASE(make_data_from_hmmer) {
//	const auto hmmer_data   = hmmer_scores_file::parse_hmmer_scores_file( hmmsearch_results_file );
//	const auto is_pos       = make_label_pair_is_positive( labelled_pair_list );
//	const auto full_results = value_lists_of_hmmer_data( hmmer_data, is_pos );
//	const auto the_set      = make_score_classn_value_results_set( full_results );
//
//	BOOST_LOG_TRIVIAL( warning ) << "The HMMER results don't cover all pairs and the code doesn't yet handle that";
//
//	auto rng = mt19937{ random_device{}() };
//	write_to_svm_light_data_files(
//		the_set,
//		hmmer_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//BOOST_AUTO_TEST_CASE(make_data_from_ssap_and_prc) {
//	const auto ssap_data         = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file     );
//	const auto prc_data          = prc_scores_file::parse_prc_scores_file_fancy   ( comp_prc_results_file );
//	const auto is_pos            = make_label_pair_is_positive( labelled_pair_list );
//	const auto plot_results_ssap = value_lists_of_ssap_data( ssap_data, is_pos );
//	const auto plot_results_prc  = value_lists_of_prc_data ( prc_data,  is_pos );
//	const auto full_results      = join_value_lists( plot_results_ssap, plot_results_prc );
//	const auto the_set           = make_score_classn_value_results_set( full_results );
//	auto rng = mt19937{ random_device{}() };
//
////	const auto scalings = get_value_list_scalings( the_set );
////	for (const auto &data_gubbins : the_set) {
////		cerr << "Series is  : " << data_gubbins.get_name() << endl;
////	}
////	for (const auto &scaling : scalings) {
////		cerr << "Scaling is : " << scaling << endl;
////	}
//
//	write_to_svm_light_data_files(
//		the_set,
//		ssap_and_prc_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//BOOST_AUTO_TEST_CASE(make_data_from_ssap_and_hmmer) {
//	const auto ssap_data          = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file      );
//	const auto hmmer_data         = hmmer_scores_file::parse_hmmer_scores_file     ( hmmsearch_results_file );
//	const auto is_pos             = make_label_pair_is_positive( labelled_pair_list );
//	const auto plot_results_ssap  = value_lists_of_ssap_data ( ssap_data,  is_pos );
//	const auto plot_results_hmmer = value_lists_of_hmmer_data( hmmer_data, is_pos );
//	const auto full_results       = join_value_lists( plot_results_ssap, plot_results_hmmer );
//	const auto the_set            = make_score_classn_value_results_set( full_results );
//	auto rng = mt19937{ random_device{}() };
//
//	write_to_svm_light_data_files(
//		the_set,
//		ssap_and_hmmer_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//BOOST_AUTO_TEST_CASE(make_data_from_ssap_prc_and_multi) {
//	const auto ssap_data         = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file        );
//	const auto prc_data          = prc_scores_file::parse_prc_scores_file_fancy   ( comp_prc_results_file    );
//	const auto multi_struc_data  = parse_multi_struc_scores                       ( multi_struc_results_file );
//
//	const auto is_pos            = make_label_pair_is_positive( labelled_pair_list );
//
//	const auto plot_results_ssap  = value_lists_of_ssap_data ( ssap_data,        is_pos );
//	const auto plot_results_prc   = value_lists_of_prc_data  ( prc_data,         is_pos );
//	const auto plot_results_multi = value_lists_of_multi_data( multi_struc_data, is_pos );
//
//	const auto full_results      = join_value_lists( plot_results_ssap, plot_results_prc, plot_results_multi );
//
//	const auto the_set           = make_score_classn_value_results_set( full_results );
//	auto rng = mt19937{ random_device{}() };
//	write_to_svm_light_data_files(
//		the_set,
//		ssap_prc_and_multi_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//BOOST_AUTO_TEST_CASE(make_data_from_ssap_prc_and_hmmer) {
//	const auto ssap_data          = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file      );
//	const auto prc_data           = prc_scores_file::parse_prc_scores_file_fancy   ( comp_prc_results_file  );
//	const auto hmmer_data         = hmmer_scores_file::parse_hmmer_scores_file     ( hmmsearch_results_file );
//	const auto is_pos             = make_label_pair_is_positive( labelled_pair_list );
//	const auto plot_results_ssap  = value_lists_of_ssap_data ( ssap_data, is_pos );
//	const auto plot_results_prc   = value_lists_of_prc_data  ( prc_data,  is_pos );
//	const auto plot_results_hmmer = value_lists_of_hmmer_data( hmmer_data, is_pos );
//	const auto full_results       = join_value_lists( plot_results_ssap, plot_results_prc, plot_results_hmmer );
//	const auto the_set            = make_score_classn_value_results_set( full_results );
//	auto rng = mt19937{ random_device{}() };
//	write_to_svm_light_data_files(
//		the_set,
//		ssap_prc_and_hmmer_svm_data_stem,
//		5,
//		rng,
//		0.5
//	);
//}

//

// BOOST_AUTO_TEST_CASE(plot_svm_from_ssap) {
// 	const auto ssap_data         = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file );
// 	const auto is_pos            = make_label_pair_is_positive( labelled_pair_list );
// 	const auto plot_results_ssap = value_lists_of_ssap_data   ( ssap_data, is_pos  );
// 	const auto plot_results_svm  = score_classn_value_list_vec{
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/ssap_rmsd.1.rbf_gamma_1_c_5.preds", "svm1" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/ssap_rmsd.2.rbf_gamma_1_c_5.preds", "svm2" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/ssap_rmsd.3.rbf_gamma_1_c_5.preds", "svm3" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/ssap_rmsd.4.rbf_gamma_1_c_5.preds", "svm4" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/ssap_rmsd.5.rbf_gamma_1_c_5.preds", "svm5" )
// 	};

// 	plot_both(
// 		join_value_lists ( plot_results_ssap,  plot_results_svm  ),
// 		( graphs_dir / "ssap_plus_svm" ).string(),
// 		join_series_lists( ssap_series_list(), svm_series_list() )
// 	);
// }

// BOOST_AUTO_TEST_CASE(plot_svm_from_prc) {
// 	const auto prc_data         = prc_scores_file::parse_prc_scores_file_fancy( comp_prc_results_file );
// 	const auto is_pos           = make_label_pair_is_positive( labelled_pair_list );
// 	const auto plot_results_prc = value_lists_of_prc_data( prc_data, is_pos );
// 	const auto plot_results_svm = score_classn_value_list_vec{
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/prc_rmsd.1.rbf_gamma_1_c_5.preds", "svm1" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/prc_rmsd.2.rbf_gamma_1_c_5.preds", "svm2" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/prc_rmsd.3.rbf_gamma_1_c_5.preds", "svm3" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/prc_rmsd.4.rbf_gamma_1_c_5.preds", "svm4" ),
// 		read_svmlight_predictions_file( "/tmp/svm_malarkey/prc_rmsd.5.rbf_gamma_1_c_5.preds", "svm5" )
// 	};

// 	plot_both(
// 		join_value_lists ( plot_results_prc,  plot_results_svm  ),
// 		( graphs_dir / "prc_plus_svm" ).string(),
// 		join_series_lists( prc_series_list(), svm_series_list() )
// 	);
// }


//BOOST_AUTO_TEST_CASE(plot_svm_from_hmmer) {
////	const auto hmmer_scan_data    = hmmer_scores_file::parse_hmmer_scores_file( hmmscan_results_file );
//	const auto hmmer_scearch_data = hmmer_scores_file::parse_hmmer_scores_file( hmmsearch_results_file );
//	const auto is_pos             = make_label_pair_is_positive( labelled_pair_list );
////	const auto plot_results_hmmer_scan   = value_lists_of_hmmer_data( hmmer_scan_data,    is_pos );
//	const auto plot_results_hmmer_search = value_lists_of_hmmer_data( hmmer_scearch_data, is_pos );
//
//	plot_both(
////		plot_results_hmmer_scan,
//		plot_results_hmmer_search,
////		( graphs_dir / "hmmer_scan" ).string(),
//		( graphs_dir / "hmmer_search" ).string(),
//		hmmer_series_list()
//	);
//}

//BOOST_AUTO_TEST_CASE(plot_svm_from_ssap_and_prc) {
//	const auto ssap_data          = ssap_scores_file::parse_ssap_scores_file_simple( ssap_results_file     );
//	const auto prc_data           = prc_scores_file::parse_prc_scores_file_fancy   ( comp_prc_results_file );
//	const auto is_pos             = make_label_pair_is_positive( labelled_pair_list );
//	const auto plot_results_ssap  = value_lists_of_ssap_data( ssap_data,           is_pos );
//	const auto plot_results_prc   = value_lists_of_prc_data ( prc_data,            is_pos );
//	const auto plot_results_magic = magic_function_data     ( ssap_data, prc_data, is_pos );
//
//	const auto plot_results_svm__prc_alone = score_classn_value_list_vec{
//		read_svmlight_predictions_file( prc_svm_data_dir                / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__prc_alone1" ),
//		read_svmlight_predictions_file( prc_svm_data_dir                / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__prc_alone2" ),
//		read_svmlight_predictions_file( prc_svm_data_dir                / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__prc_alone3" ),
//		read_svmlight_predictions_file( prc_svm_data_dir                / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__prc_alone4" ),
//		read_svmlight_predictions_file( prc_svm_data_dir                / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__prc_alone5" )
//	};
//
//	const auto plot_results_svm__ssap_alone = score_classn_value_list_vec{
//		read_svmlight_predictions_file( ssap_svm_data_dir               / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__ssap_alone1" ),
//		read_svmlight_predictions_file( ssap_svm_data_dir               / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__ssap_alone2" ),
//		read_svmlight_predictions_file( ssap_svm_data_dir               / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__ssap_alone3" ),
//		read_svmlight_predictions_file( ssap_svm_data_dir               / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__ssap_alone4" ),
//		read_svmlight_predictions_file( ssap_svm_data_dir               / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__ssap_alone5" )
//	};
//
//	const auto plot_results_svm__hmmer_alone_prc = score_classn_value_list_vec{
//		read_svmlight_predictions_file( hmmer_svm_data_dir              / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__hmmer_alone1" ),
//		read_svmlight_predictions_file( hmmer_svm_data_dir              / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__hmmer_alone2" ),
//		read_svmlight_predictions_file( hmmer_svm_data_dir              / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__hmmer_alone3" ),
//		read_svmlight_predictions_file( hmmer_svm_data_dir              / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__hmmer_alone4" ),
//		read_svmlight_predictions_file( hmmer_svm_data_dir              / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__hmmer_alone5" )
//	};
//
//	const auto plot_results_svm__ssap_prc   = score_classn_value_list_vec{
//		read_svmlight_predictions_file( ssap_and_prc_svm_data_dir       / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__ssap_prc1" ),
//		read_svmlight_predictions_file( ssap_and_prc_svm_data_dir       / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__ssap_prc2" ),
//		read_svmlight_predictions_file( ssap_and_prc_svm_data_dir       / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__ssap_prc3" ),
//		read_svmlight_predictions_file( ssap_and_prc_svm_data_dir       / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__ssap_prc4" ),
//		read_svmlight_predictions_file( ssap_and_prc_svm_data_dir       / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__ssap_prc5" )
//	};
//
//	const auto plot_results_svm__ssap_hmmer   = score_classn_value_list_vec{
//		read_svmlight_predictions_file( ssap_and_hmmer_svm_data_dir     / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__ssap_hmmer1" ),
//		read_svmlight_predictions_file( ssap_and_hmmer_svm_data_dir     / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__ssap_hmmer2" ),
//		read_svmlight_predictions_file( ssap_and_hmmer_svm_data_dir     / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__ssap_hmmer3" ),
//		read_svmlight_predictions_file( ssap_and_hmmer_svm_data_dir     / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__ssap_hmmer4" ),
//		read_svmlight_predictions_file( ssap_and_hmmer_svm_data_dir     / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__ssap_hmmer5" )
//	};
//
//	const auto plot_results_svm__ssap_prc_hmmer = score_classn_value_list_vec{
//		read_svmlight_predictions_file( ssap_prc_and_hmmer_svm_data_dir / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__ssap_prc_hmmer1" ),
//		read_svmlight_predictions_file( ssap_prc_and_hmmer_svm_data_dir / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__ssap_prc_hmmer2" ),
//		read_svmlight_predictions_file( ssap_prc_and_hmmer_svm_data_dir / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__ssap_prc_hmmer3" ),
//		read_svmlight_predictions_file( ssap_prc_and_hmmer_svm_data_dir / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__ssap_prc_hmmer4" ),
//		read_svmlight_predictions_file( ssap_prc_and_hmmer_svm_data_dir / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__ssap_prc_hmmer5" )
//	};
//
//	const auto plot_results_svm__ssap_prc_multi = score_classn_value_list_vec{
//		read_svmlight_predictions_file( ssap_prc_and_multi_svm_data_dir / "svm_data.1.rbf_gamma_1_c_5.preds", "svm__ssap_prc_multi1" ),
//		read_svmlight_predictions_file( ssap_prc_and_multi_svm_data_dir / "svm_data.2.rbf_gamma_1_c_5.preds", "svm__ssap_prc_multi2" ),
//		read_svmlight_predictions_file( ssap_prc_and_multi_svm_data_dir / "svm_data.3.rbf_gamma_1_c_5.preds", "svm__ssap_prc_multi3" ),
//		read_svmlight_predictions_file( ssap_prc_and_multi_svm_data_dir / "svm_data.4.rbf_gamma_1_c_5.preds", "svm__ssap_prc_multi4" ),
//		read_svmlight_predictions_file( ssap_prc_and_multi_svm_data_dir / "svm_data.5.rbf_gamma_1_c_5.preds", "svm__ssap_prc_multi5" )
//	};
//
//	const auto plot_results_svm__full = score_classn_value_list_vec{
//		read_svmlight_predictions_file( ssap_and_prc_svm_data_dir       / "svm_data.full.rbf_gamma_1_c_5.preds", "svm__full" )
//	};
//
//	plot_both(
//		join_value_lists ( plot_results_ssap,  plot_results_prc,  plot_results_magic,  plot_results_svm__hmmer_alone_prc, plot_results_svm__prc_alone,  plot_results_svm__ssap_alone,  plot_results_svm__ssap_prc_hmmer, plot_results_svm__ssap_prc_multi,  plot_results_svm__ssap_hmmer,  plot_results_svm__ssap_prc,  plot_results_svm__full  ),
//		( graphs_dir / "ssap_and_prc_plus_svm" ).string(),
//		join_series_lists( ssap_series_list(), prc_series_list(), magic_series_list(), svm__hmmer_alone_series_list(),    svm__prc_alone_series_list(), svm__ssap_alone_series_list(), svm__ssap_prc_hmmer_series_list(),svm__ssap_prc_multi_series_list(), svm__ssap_hmmer_series_list(), svm__ssap_prc_series_list(), svm__full_series_list() )
//	);
//}

// cd /tmp
// wget "http://download.joachims.org/svm_light/current/svm_light.tar.gz"
// mkdir svm_light
// tar --directory svm_light -zxvf svm_light.tar.gz
// cd svm_light
// chmod -x *
// make -j 3
//
// # Default: -t 2 -g 1 -c 5
// # Explore: -g varying over  {  2^-20,  2^-19,  2^-18,... ,  2^5 }
// # Explore: -c varying over  { 10^-20, 10^-19, 10^-18,... , 10^5 }
//
//
// seq 5 | xargs -P 4 -I VAR ./svm_learn -t 2 -g 1 -c 5       /tmp/svm_malarkey/ssap_rmsd.VAR.train /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.model
// seq 5 | xargs -P 4 -I VAR ./svm_classify                   /tmp/svm_malarkey/ssap_rmsd.VAR.test  /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.model                                                      /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.preds.raw
// seq 5 | xargs      -I VAR /bin/tcsh -c "paste -d ' '       /tmp/svm_malarkey/ssap_rmsd.VAR.test  /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.preds.raw > /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.preds.raw.paste"
// echo '{print $(NF-2)" "$(NF-1)" "$(NF)" "$1}' > /tmp/ap
// seq 5 | xargs      -I VAR /bin/tcsh -c "awk --file /tmp/ap /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.preds.raw.paste > /tmp/svm_malarkey/ssap_rmsd.VAR.rbf_gamma_1_c_5.preds"

BOOST_AUTO_TEST_SUITE_END()





//BOOST_AUTO_TEST_SUITE(score_classn_value_results_set_test_suite)
//
///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(ssap) {
//	auto som2_results = make_score_classn_value_list_vec( get_som2_scores() );
//	plot_both(
//		som2_results,
//		"/export/people/ucbctnl/sab_graphs/1.ssap",
//		{ { "ssap", str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } } }
//	);
//
////	plot_roc(
////		classn_stat_plotter(),
////		"/export/people/ucbctnl/sab_graphs/1.ssap.roc",
////		som2_results,
////		make_standard_score_roc_plotter_spec( { { "ssap", str_opt{ R"( linetype 1 linecolor rgb "#000000" )" } } } )
////	);
////
////	plot_precision_recall(
////		classn_stat_plotter(),
////		"/export/people/ucbctnl/sab_graphs/1.ssap.precision_recall",
////		som2_results,
////		make_standard_score_precision_recall_plotter_spec( { { "ssap", str_opt{ R"( linetype 1 linecolor rgb "#000000" )" } } } )
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
//			{ "ssap",                                                 str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "shorter_protein_length",                               str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "geometric_mean_protein_length",                        str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length", str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                  str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "longer_protein_length",                                str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } }
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
//			{ "ssap",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "GSAS",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "SAS",                                                                                    str_opt{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "SIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "num_aligned_residues",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    str_opt{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
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
//			{ "ssap",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
////			{ "GSAS",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
////			{ "SAS",                                                                                    str_opt{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#2222ff" linewidth 3 )" } },
////			{ "MIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
////			{ "SIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
////			{ "sequence_id.identity",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
////			{ "lDDT.threshold[STD_MEAN]",                                                               str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
////			{ "num_aligned_residues",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
////			{ "overlap.num_aligned_residues.longer_protein_length",                                     str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
////			{ "dRMSD",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
////			{ "shorter_protein_length",                                                                 str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
////			{ "SI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
////			{ "MI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
////			{ "RMSD",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
////			{ "geometric_mean_protein_length",                                                          str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
////			{ "overlap.shorter_protein_length.longer_protein_length",                                   str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
////			{ "mean_protein_length",                                                                    str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
////			{ "sequence_id.blosum62",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
////			{ "longer_protein_length",                                                                  str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
////			{ "overlap.num_aligned_residues.shorter_protein_length",                                    str_opt{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//			{ "Native TM-align",                                                                        str_opt{ R"( linetype 1 linecolor rgb "#000088" linewidth 3 )" } }
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
//			{ "ssap",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
////			{ "GSAS",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
////			{ "SAS",                                                                                    str_opt{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#2222ff" linewidth 3 )" } },
////			{ "MIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
////			{ "SIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
////			{ "sequence_id.identity",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
////			{ "lDDT.threshold[STD_MEAN]",                                                               str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
////			{ "num_aligned_residues",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
////			{ "overlap.num_aligned_residues.longer_protein_length",                                     str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
////			{ "dRMSD",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
////			{ "shorter_protein_length",                                                                 str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
////			{ "SI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
////			{ "MI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
////			{ "RMSD",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
////			{ "geometric_mean_protein_length",                                                          str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
////			{ "overlap.shorter_protein_length.longer_protein_length",                                   str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
////			{ "mean_protein_length",                                                                    str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
////			{ "sequence_id.blosum62",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
////			{ "longer_protein_length",                                                                  str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
////			{ "overlap.num_aligned_residues.shorter_protein_length",                                    str_opt{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
//			{ "Native TM-align",                                                                        str_opt{ R"( linetype 1 linecolor rgb "#000088" linewidth 3 )" } }
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
//			{ "ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified", str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70].cb_atoms",                                            str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "ssap",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "SAS.select_best_score_percent[70].cb_atoms",                                             str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70].cb_atoms",                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70]",                                                 str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70]",                                                     str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_best_score_percent[70]",                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.cb_atoms",                                                                          str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01].cb_atoms",                                                   str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01].cb_atoms",                                               str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01]",                                                            str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "structal.cb_atoms",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.select_min_score[0.01].cb_atoms",                                               str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.cb_atoms",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01].cb_atoms",                                                    str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01]",                                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.cb_atoms",                                                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01]",                                                             str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS",                                                                                    str_opt{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score.select_min_score[0.01]",                                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX.select_best_score_percent[70].cb_atoms",                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX.select_best_score_percent[70].cb_atoms",                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "MIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "lDDT.threshold[4]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[2]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "lDDT.threshold[1]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "num_aligned_residues",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "lDDT.threshold[0.5]",                                                                    str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    str_opt{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
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
//			{ "SVM01",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM02",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM03",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM04",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM05",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM06",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM07",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM08",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM09",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM10",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//
//
//			{ "ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified", str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70].cb_atoms",                                            str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "ssap",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "SAS.select_best_score_percent[70].cb_atoms",                                             str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70].cb_atoms",                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70]",                                                 str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70]",                                                     str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_best_score_percent[70]",                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.cb_atoms",                                                                          str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01].cb_atoms",                                                   str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01].cb_atoms",                                               str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01]",                                                            str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "structal.cb_atoms",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.select_min_score[0.01].cb_atoms",                                               str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.cb_atoms",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01].cb_atoms",                                                    str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01]",                                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.cb_atoms",                                                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01]",                                                             str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS",                                                                                    str_opt{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score.select_min_score[0.01]",                                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX.select_best_score_percent[70].cb_atoms",                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX.select_best_score_percent[70].cb_atoms",                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "MIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "lDDT.threshold[4]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[2]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "lDDT.threshold[1]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "num_aligned_residues",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "lDDT.threshold[0.5]",                                                                    str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    str_opt{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
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
//			{ "SVM01",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM02",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM03",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM04",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM05",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM06",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM07",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM08",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM09",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//			{ "SVM10",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#FFCC00" notitle     )" } },
//
//
//			{ "ssap.cb_atoms.high_accuracy.num_excluded_on_sides:10.distance_score_formula_simplified", str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70].cb_atoms",                                            str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "ssap",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "black"   linewidth 3 )" } },
//			{ "SAS.select_best_score_percent[70].cb_atoms",                                             str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70].cb_atoms",                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_best_score_percent[70]",                                                 str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_best_score_percent[70]",                                                     str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_best_score_percent[70]",                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.cb_atoms",                                                                          str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01].cb_atoms",                                                   str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01].cb_atoms",                                               str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS.select_min_score[0.01]",                                                            str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "GSAS",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#00ff00"             )" } },
//			{ "structal.cb_atoms",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.select_min_score[0.01].cb_atoms",                                               str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score.cb_atoms",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01].cb_atoms",                                                    str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "structal.select_min_score[0.01]",                                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.cb_atoms",                                                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS.select_min_score[0.01]",                                                             str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SAS",                                                                                    str_opt{ R"( linetype 1 linecolor rgb "#0000ff"             )" } },
//			{ "structal",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#ff0000"             )" } },
//			{ "TM-score.select_min_score[0.01]",                                                        str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "TM-score",                                                                               str_opt{ R"( linetype 1 linecolor rgb "#01fffe"             )" } },
//			{ "MIMAX.select_best_score_percent[70].cb_atoms",                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX.select_best_score_percent[70].cb_atoms",                                           str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "MIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffa6fe"             )" } },
//			{ "lDDT.threshold[4]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "SIMAX",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ffdb66"             )" } },
//			{ "sequence_id.identity",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#006401"             )" } },
//			{ "lDDT.threshold[2]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "lDDT.threshold[STD_MEAN]",                                                               str_opt{ R"( linetype 1 linecolor rgb "#010067"             )" } },
//			{ "lDDT.threshold[1]",                                                                      str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "num_aligned_residues",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#95003a"             )" } },
//			{ "lDDT.threshold[0.5]",                                                                    str_opt{ R"( linetype 2 linecolor rgb "#cccccc"             )" } },
//			{ "overlap.num_aligned_residues.longer_protein_length",                                     str_opt{ R"( linetype 1 linecolor rgb "#007db5"             )" } },
//			{ "dRMSD",                                                                                  str_opt{ R"( linetype 1 linecolor rgb "#ff00f6"             )" } },
//			{ "shorter_protein_length",                                                                 str_opt{ R"( linetype 1 linecolor rgb "#774d00"             )" } },
//			{ "SI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#90fb92"             )" } },
//			{ "MI",                                                                                     str_opt{ R"( linetype 1 linecolor rgb "#0076ff"             )" } },
//			{ "RMSD",                                                                                   str_opt{ R"( linetype 1 linecolor rgb "#d5ff00"             )" } },
//			{ "geometric_mean_protein_length",                                                          str_opt{ R"( linetype 1 linecolor rgb "#ff937e"             )" } },
//			{ "overlap.shorter_protein_length.longer_protein_length",                                   str_opt{ R"( linetype 1 linecolor rgb "#6a826c"             )" } },
//			{ "mean_protein_length",                                                                    str_opt{ R"( linetype 1 linecolor rgb "#ff029d"             )" } },
//			{ "sequence_id.blosum62",                                                                   str_opt{ R"( linetype 1 linecolor rgb "#fe8900"             )" } },
//			{ "longer_protein_length",                                                                  str_opt{ R"( linetype 1 linecolor rgb "#7a4782"             )" } },
//			{ "overlap.num_aligned_residues.shorter_protein_length",                                    str_opt{ R"( linetype 1 linecolor rgb "#7e2dd2"             )" } }
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
//			{ "ssap",                 str_opt{ R"( linetype 1 linecolor rgb "red"    linewidth 3 )" } },
//			{ "dssp_based_cath-ssap", str_opt{ R"( linetype 1 linecolor rgb "blue"   linewidth 3 )" } },
//			{ "SSAP_r17937",          str_opt{ R"( linetype 1 linecolor rgb "green"  linewidth 3 )" } },
//			{ "SSAP_r7906",           str_opt{ R"( linetype 1 linecolor rgb "orange" linewidth 3 )" } }
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
////			using bool_str_pair_opt = optional<pair<bool, string>>;
////			if ( file_size( arg_file ) < 100 || icontains( arg_file.string(), "stderr" ) ) {
////				return bool_str_pair_opt( none );
////			}
////
////			const auto file_stem = arg_file.stem().string();
////			return bool_str_pair_opt{ make_pair(
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

