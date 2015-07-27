/// \file
/// \brief The all_vs_all test suite

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London.
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

#include <boost/log/trivial.hpp> // ***** TEMPORARY *****
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/sort.hpp> // ***** TEMPORARY *****
#include <boost/range/irange.hpp>

#include "common/boost_addenda/range/adaptor/limited.h" // ***** TEMPORARY *****
#include "common/chrono/duration_to_seconds_string.h"
#include "common/file/simple_file_read_write.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "scan/scan_action/record_scores_scan_action.h"
#include "scan/scan_tools/all_vs_all.h"
#include "scan/scan_tools/scan_metrics.h"
#include "score/pair_scatter_plotter/pair_scatter_plotter.h"  // ***** TEMPORARY *****
#include "structure/protein/protein_source_file_set/protein_source_from_pdb.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"
#include "test/global_test_constants.h"

#include <chrono>

using namespace cath;
using namespace cath::common;
using namespace cath::scan;
using namespace cath::score;
using namespace std;
using namespace std::chrono;

using boost::irange;
using boost::numeric_cast;
using boost::range::sort;

namespace cath {
	namespace test {

		/// \brief The all_vs_all_test_suite_fixture to assist in testing all_vs_all
		struct all_vs_all_test_suite_fixture : protected global_test_constants {
		protected:
			~all_vs_all_test_suite_fixture() noexcept = default;



			str_str_pair_doub_map read_ssap_scores() const {
				const auto raw_scores = read_file<str_str_pair_doub_pair>( "ssap_score_summary_file.txt" );
				str_str_pair_doub_map ssap_scores;
				for (const auto &x : raw_scores) {
					ssap_scores.insert( x );
				}
				return ssap_scores;
			}

			/// \brief TODOCUMENT
			double calc_plot_colour_component(const size_size_pair &arg_length_range,
			                                  const size_t         &arg_length
			                                  ) {
				const auto &min = arg_length_range.first;
				const auto &max = arg_length_range.second;
				if ( min > max || arg_length < min || arg_length > max ) {
					BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calc_plot_colour_component() for invalid length/range"));
				}

				const auto length_less_min = numeric_cast<double>( arg_length - min );
				const auto max_less_min    = numeric_cast<double>(        max - min );
				return ( max == min) ? 0.0
				                     : ( length_less_min / max_less_min );
			}
		};

	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(all_vs_all_test_suite, cath::test::all_vs_all_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	const auto full_ssap_scores = read_ssap_scores();

	ostringstream parse_ss;
//	const auto ids             = str_vec{ "1a04A02", "1c0pA01", "1fseB00", "1h3fA01", "1hdoA00", "1n3lA01", "1r6xA02" };
//	const auto ids             = str_vec{ "156b000", "1a02F00", "1a04A01", "1a04A02", "1a0cA00", "1a0iA02", "1a0iA03", "1a1hA01", "1a2oA01", "1a2pA00", "1a2xB00", "1a32A00", "1a34A00", "1a3qA01", "1a4iB01", "1a5tA01", "1a62A02", "1a6bB00", "1a6cA01", "1a76A01", "1a7jA00", "1a7sA02", "1a8dA01", "1a9xA04", "1a9xA08", "1a9xB02", "1aa0A00", "1ac0A00", "1ac5A00", "1aepA00", "1afpA00", "1aooA00", "1arbA01", "1atgA01", "1au7A02", "1avcA07", "1avwB00", "1avyA00", "1avyB00", "1b06A01", "1b25A01", "1b4pA01", "1b8xA03", "1b9mB02", "1bboA02", "1bd8A00", "1be3I00", "1bg5A03", "1bx7A00", "1cf7B00", "1cl4A00", "1cptA00", "1cw5A00", "1d0cA01", "1devB00", "1dk5B01", "1f8vE00", "1fcyA00", "1fmjA00", "1fseB00", "1g3jD00", "1g4mA01", "1g8kA02", "1giqA01", "1go3F02", "1gpzA02", "1h4aX01", "1h59B00", "1h6gA02", "1h6wA02", "1hfeL02", "1hp1A02", "1hwmA02", "1hykA00", "1i4oC00", "1i7wD00", "1ifcA00", "1j0hA02", "1kk1A01", "1kyiS02", "1kyqA02", "1m2tA02", "1m32A02", "1m7vA01", "1m9oA00", "1ma1B01", "1mhwC00", "1mo9A01", "1mslA02", "1mtyB00", "1my5A00", "1my7A00", "1n7dA02", "1ncqD00", "1nn6A01", "1od3A00", "1ojjA00", "1okiA01", "1p3cA02", "1p9oA00", "1pbgA00", "1pc3A01", "1pzlA00", "1qfoC00", "1qo0D02", "1r6vA03", "1rqgA04", "1rr7A02", "1rrkA02", "1ryp100", "1s70B01", "1sc6A01", "1si5H01", "1sl6A00", "1sqjB02", "1szbA02", "1t6lA00", "1tnsA00", "1u7lA01", "1ufmA00", "1ukcB00", "1uqtA01", "1uwwA00", "1vsgA01", "1vyrA00", "1w6gA03", "1wfqA01", "1wkqB00", "1x6mA00", "1xksA00", "1y8qB01", "1yarO01", "1yewA01", "1ytqA01", "1zodA02", "2a4dA00", "2ayzA00", "2bs2C00", "2c4jA01", "2cz9A02", "2dlkA02", "2e6iA00", "2fe5A00", "2fmpA04", "2g3aA01", "2grjA00", "2gy5A03", "2h7cB00", "2i7dA01", "2ieeA01", "2ij2B00", "2iy5A00", "2j7jA03", "2jn4A00", "2mev400", "2nv0A00", "2o01N01", "2o2cC01", "2o55A00", "2ob5A00", "2odkA00", "2oq1A02", "2ph1A00", "2pw9A02", "2pw9C02", "2qcvA01", "2qjpB02", "2qjyB02", "2r1fA01", "2rh0A01", "2uzyB02", "2v3iA00", "2vanA03", "2vglB00", "2vpzC00", "2w15A00", "2wjwA01", "2xblD00", "2xw9A02", "2y1eA01", "3arcO02", "3bixA00", "3bmvA01", "3broD00", "3c2uA01", "3cb2B01", "3cjwA00", "3d64A01", "3f94A02", "3fuyA00", "3g02B00", "3gixA00", "3i07B01", "3iauB01", "3im9A01", "3kvnA01", "3laeA00", "3ljkA01", "3nwnA00", "3o2zH00", "3oaeA00", "3p4tA03", "3qzbA00" };
//	const auto ids             = str_vec{ "1a04A01", "1a2pA00", "1a2xB00", "1a34A00", "1a7jA00", "1aooA00", "1avyA00", "1b25A01", "1b4pA01", "1be3I00", "1d0cA01", "1devB00", "1f8vE00", "1fcyA00", "1g3jD00", "1g8kA02", "1hfeL02", "1i4oC00", "1i7wD00", "1kyiS02", "1ncqD00", "1ojjA00", "1p9oA00", "1rqgA04", "1sqjB02", "1u7lA01", "1vsgA01", "1vyrA00", "1w6gA03", "1wkqB00", "1x6mA00", "1y8qB01", "1zodA02", "2c4jA01", "2gy5A03", "2ieeA01", "2jn4A00", "2nv0A00", "2o01N01", "2odkA00", "2ph1A00", "2pw9A02", "2qcvA01", "2qjpB02", "2qjyB02", "2rh0A01", "2uzyB02", "2vglB00", "2vpzC00", "2xw9A02", "3c2uA01", "3cjwA00", "3f94A02", "3g02B00", "3iauB01", "3laeA00", "3nwnA00", "3o2zH00", "1a04A02", "1fseB00" };
	const auto ids             = str_vec{ "1my7A00", "1my5A00", "2qjyB02", "2qjpB02", "2pw9A02", "2pw9C02", "2c4jA01", "1b4pA01", "2fmpA04", "2vanA03", "1okiA01", "1ytqA01", "1b06A01", "1ma1B01", "1a7sA02", "2xw9A02", "1avyB00", "1avyA00", "1m2tA02", "1hwmA02", "1d0cA01", "1m7vA01", "1a1hA01", "2j7jA03", "1a04A02", "1fseB00", "1fcyA00", "1pzlA00", "1avcA07", "1dk5B01", "1bd8A00", "1s70B01", "1atgA01", "1pc3A01", "1a2oA01", "2ayzA00", "1au7A02", "1rr7A02", "1arbA01", "1si5H01", "1ufmA00", "1a9xB02", "2nv0A00", "1aepA00", "1h6gA02", "1a4iB01", "1sc6A01", "2y1eA01", "1cf7B00", "1a32A00", "1go3F02", "3broD00", "1tnsA00", "2xblD00", "1a3qA01", "1g4mA01", "1a04A01", "2wjwA01", "1a02F00", "1mslA02" };
//	const auto load_files_starttime = high_resolution_clock::now();
	const auto proteins             = read_proteins_from_files( protein_source_from_pdb(), TEST_SOURCE_DATA_DIR(), ids, parse_ss );
//	const auto load_files_duration  = high_resolution_clock::now() - load_files_starttime;

	const auto query_length_min_max = min_max_protein_length( proteins );
	const auto match_length_min_max = min_max_protein_length( proteins );


	const auto do_stuff_start = high_resolution_clock::now();
	const auto the_scan_action_and_metrics = all_vs_all{}.perform_scan( proteins, proteins );
	const auto do_stuff_durn  = high_resolution_clock::now() - do_stuff_start;

//	the_scan_action_and_metrics

	BOOST_LOG_TRIVIAL( warning ) << "Did stuff - took " << durn_to_seconds_string        ( do_stuff_durn )
	                             << " ("                << durn_to_rate_per_second_string( do_stuff_durn ) << ")";
	const auto prot_idx_range  = irange( 0_z, proteins.size() );
	cerr << "       ";
	for (const size_t match_index : prot_idx_range) {
		cerr << "  " << ids[ match_index ];
	}
	cerr << "\n";
	vector<tuple<double, double, double, double, double>> results_for_plotting;
//	doub_doub_pair_vec results_for_plotting;
//	size_size_doub_tpl_vec query_match_score_entries;
	for (const size_t query_index : prot_idx_range) {
		const auto &query_id = ids[ query_index ];
		cerr << query_id;
		for (const size_t match_index : prot_idx_range) {
			const auto   &match_id         = ids[ match_index ];
			const size_t  longest_length   = max( proteins[ query_index ].get_length(), proteins[ match_index ].get_length() );
			const double  raw_score        = the_scan_action_and_metrics.first.get_score( query_index, match_index );
			const double  normalised_score = raw_score / numeric_cast<double>( ( longest_length - 10 ) * ( longest_length - 11 ) );
			cerr << "WARNING: ******** ASSUMING THAT THE MINIMUM_INDEX_DISTANCE IS 11 ***********\n";
			const double  logged_score     = 1.0 + log( normalised_score ) / log( 50000.0 );
			const double  final_score      = ( raw_score != 0.0 ) ? 100.0 * logged_score
			                                                     :  -15.0;
			cerr << " " << setw( 8 ) << right << final_score;
			const auto &ssap_score = full_ssap_scores.at( make_pair( ids[ query_index ], ids[ match_index ] ) );
//			query_match_score_entries.emplace_back( query_index, match_index, final_score );
			results_for_plotting.emplace_back(
				ssap_score,
				final_score,
				calc_plot_colour_component( query_length_min_max, proteins[ query_index ].get_length() ),
				calc_plot_colour_component( match_length_min_max, proteins[ match_index ].get_length() ),
				0.0
			);
			if ( ssap_score  > 75.0 && final_score < 25.0 ) {
				cerr << "TO INVESTIGATE FURTHER : " << query_id << " " << match_id << " SSAP: " << ssap_score << " SNAP: " << final_score << "\n";
			}
		}
		cerr << "\n";
	}

//	rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)
//	xrgb(r,g,b) = (g-b)/255. * cos(30.)
//	yrgb(r,g,b) = r/255. - (g+b)/255. * sin(30.)
//	GPFUN_rgb = "rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)"
//	GPFUN_xrgb = "xrgb(r,g,b) = (g-b)/255. * cos(30.)"
//	GPFUN_yrgb = "yrgb(r,g,b) = r/255. - (g+b)/255. * sin(30.)"
//	splot 'rgb_variable.dat' using 1:2:3:(rgb($1,$2,$3)) with points pt 7 ps 4 lc rgb variable,       '' using 1:2:3:(sprintf("0x%x",rgb($1,$2,$3))) with labels left offset 1 notitle

//rgb(r,g,b) = int(255*(1-r))*65536 + int(255*(1-g))*256 + int(255*(1 - 0.5 * r - 0.5 * g))
//set xtics -10,10
//set ytics -10,10
//set xrange [-15:110]
//set yrange [-15:110]
//set object 1 rectangle from -15,-15 to 110,110 fillcolor rgb "black" behind
//plot x with lines lc rgb "#B8B8B8" notitle, '/tmp/ssap_snap.data.txt' using 1:2:(rgb($3,$4,$5)) with points pointtype 7 pointsize 0.75 lc rgb variable notitle

	// set xtics -10,10
	// set ytics -10,10
	// set xrange [-15:110]
	// set yrange [-15:110]
	// plot x with lines lc rgb "#B8B8B8" notitle, '/tmp/ssap_snap.data.txt'  with points pointtype 7 pointsize 0.25 lc rgb '#000088' notitle
	pair_scatter_plotter().plot( "/tmp/ssap_snap", results_for_plotting, "SSAP", "SNAP" );

//	sort(
//		query_match_score_entries,
//		[] (const size_size_doub_tpl &x, const size_size_doub_tpl &y) {
//			return get<2>( x ) < get<2>( y );
//		}
//	);
//
//	for (const auto &x : query_match_score_entries) {
//		const auto &query_index = get<0>( x );
//		const auto &match_index = get<1>( x );
//		const auto &score       = get<2>( x );
//		if ( query_index != match_index ) {
//			const auto &query_id   = ids     [ query_index ];
//			const auto &match_id   = ids     [ match_index ];
//			const auto &query_prot = proteins[ query_index ];
//			const auto &match_prot = proteins[ match_index ];
//			cerr << query_id
//			     << " (" << right << setw(5) << query_prot.get_length() << ") "
//			     << match_id
//				 << " (" << right << setw(5) << match_prot.get_length() << ") -> "
//				 << score
//				 << "\n";
//		}
//	}

	// 1a04A02  1fseB00   80   70  87.49   67   83   30   4.77
	// 1n3lA01  1h3fA01  209  195  84.91  186   88   26   2.68
	BOOST_CHECK_EQUAL(0, 0);
}

BOOST_AUTO_TEST_SUITE_END()
