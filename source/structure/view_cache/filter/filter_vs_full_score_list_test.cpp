/// \file
/// \brief The filter_vs_full_score_list test suite

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

#include "common/file/simple_file_read_write.hpp"
#include "common/size_t_literal.hpp"
#include "score/true_pos_false_neg/classn_num_stat.hpp"
#include "score/true_pos_false_neg/classn_rate_stat.hpp"
#include "score/true_pos_false_neg/true_false_pos_neg.hpp"
#include "structure/view_cache/filter/filter_vs_full_score.hpp"
#include "structure/view_cache/filter/filter_vs_full_score_list.hpp"

#include <random>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace cath::index::filter;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The filter_vs_full_score_list_test_suite_fixture to assist in testing filter_vs_full_score_list
		struct filter_vs_full_score_list_test_suite_fixture {
		protected:
			~filter_vs_full_score_list_test_suite_fixture() noexcept = default;

		public:
		    filter_vs_full_score_list generate_random_filter_vs_full_score_list(const size_t &) const;
		};

	}
}

/// \brief Generate some data with random noise that can be used in other tests
///
/// The generated data is based on a curve `full_score = sqrt( filter_score )` over [0,1) but it's done by:
///  * randomly (uniformly) generating filter_scores in [0,1) and then squaring to get a filter score
///  * randomly (uniformly) generating some noise to add to the filter score (noise_width)
///  * before adding the noise, scaling the filter score down so that range with the noise will be [0,1)
filter_vs_full_score_list cath::test::filter_vs_full_score_list_test_suite_fixture::generate_random_filter_vs_full_score_list(const size_t &arg_num_entries ///< The number of entries the randomly generated filter_vs_full_score_list should contain
                                                                                                                              ) const {
	// Specify the width of the random noise, and calculate the scaling that it implies
	const double noise_width       = 0.1;
	const double pre_noise_scaling = 1.0 - noise_width;

	// Prepare relevant distributions and a random number generator
	uniform_real_distribution<double> full_score_distn( 0.0, 1.0 );
	uniform_real_distribution<double> noise_distn     ( -0.5 * noise_width, 0.5 * noise_width );
	default_random_engine rng{ random_device{}() };

	// Loop over the number of entries, building a vector of filter_vs_full_score objects
	filter_vs_full_score_vec new_filter_vs_full_scores;
	for (size_t entry_ctr = 0; entry_ctr < arg_num_entries; ++entry_ctr) {
		const double full_score         = full_score_distn( rng );
		const double raw_filter_score   = full_score * full_score;
		const double clean_filter_score = raw_filter_score * pre_noise_scaling + noise_width / 2.0;
		const double filter_score       = clean_filter_score + noise_distn( rng );
		new_filter_vs_full_scores.push_back( filter_vs_full_score( 100.0 * filter_score, 100.0 * full_score ) );
	}

	// Return a filter_vs_full_score_list from the vector of filter_vs_full_score objects
	return filter_vs_full_score_list( new_filter_vs_full_scores );
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(filter_vs_full_score_list_test_suite, cath::test::filter_vs_full_score_list_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
//	const double test_data_cutoff = 0.98;
//	const filter_vs_full_score_list data    = generate_random_filter_vs_full_score_list( 5000 );
//	const true_false_pos_neg        result  = filter_result_full_score_with_sensitivity ( data, 0.5, test_data_cutoff );
//	const filter_vs_full_score_vec  filters = filter_attempts_with_sensitivity( data, test_data_cutoff );
//	gnuplot_data                   ( data, path( "filter_vs_full_score_list"            ), filters );
//	gnuplot_classsn_stat_for_recall( data, path( "filter_precision_for_99_recall"       ), precision(),           test_data_cutoff );
//	gnuplot_classsn_stat_for_recall( data, path( "filter_false_positives_for_99_recall" ), false_positive_stat(), test_data_cutoff );
//
//
//	const double random_sensitivity_cutoff = 0.98;
//	const filter_vs_full_score_vec  raw_random_data = read_file<filter_vs_full_score>(path("random_filter_data.txt"));
//	const filter_vs_full_score_list random_data( raw_random_data );
//	const filter_vs_full_score_vec  random_filters = filter_attempts_with_sensitivity( random_data, random_sensitivity_cutoff );
//
//	gnuplot_data                   ( random_data, path( "random_vs_ssap_score_list"                   ), random_filters );
//	gnuplot_classsn_stat_for_recall( random_data, path( "random_filter_precision_for_99_recall"       ), precision(),           random_sensitivity_cutoff );
//	gnuplot_classsn_stat_for_recall( random_data, path( "random_filter_false_positives_for_99_recall" ), false_positive_stat(), random_sensitivity_cutoff );
//
//	const double grath_sensitivity_cutoff = 0.98;
//	const filter_vs_full_score_vec  raw_grath_data = read_file<filter_vs_full_score>(path("grath_filter_data.txt"));
//	const filter_vs_full_score_list grath_data( raw_grath_data );
//	const filter_vs_full_score_vec  grath_filters = filter_attempts_with_sensitivity( grath_data, grath_sensitivity_cutoff );
//
//	gnuplot_data                   ( grath_data, path( "grath_vs_ssap_score_list"                   ), grath_filters );
//	gnuplot_classsn_stat_for_recall( grath_data, path( "grath_filter_precision_for_99_recall"       ), precision(),           grath_sensitivity_cutoff );
//	gnuplot_classsn_stat_for_recall( grath_data, path( "grath_filter_false_positives_for_99_recall" ), false_positive_stat(), grath_sensitivity_cutoff );


	BOOST_CHECK( true );
}

//	cerr << "result is   : " << result << endl;
//	cerr << "sensitivity : " << rational_cast<double>( sensitivity().calculate( result ) ) << endl;
//	cerr << "specificity : " << rational_cast<double>( specificity().calculate( result ) ) << endl;
//	cerr << "precision   : " << rational_cast<double>( precision  ().calculate( result ) ) << endl;

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(random_generation) {
	BOOST_CHECK_EQUAL( generate_random_filter_vs_full_score_list( 30 ).size(), 30_z );
}

BOOST_AUTO_TEST_SUITE_END()
