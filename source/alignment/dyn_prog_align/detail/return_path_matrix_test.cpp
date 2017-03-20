/// \file
/// \brief The return_path_matrix test suite

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

#include "alignment/dyn_prog_align/detail/return_path_matrix.hpp"
#include "test/global_test_constants.hpp"

#include <random>

using namespace cath::align::detail;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The return_path_matrix_test_suite_fixture to assist in testing return_path_matrix
		struct return_path_matrix_test_suite_fixture {
		protected:
			~return_path_matrix_test_suite_fixture() noexcept = default;

			path_step make_random_step() const;

			/// \brief TODOCUMENT
			return_path_matrix make_random_return_path_matrix(const size_t &,
			                                                  const size_t &,
			                                                  const size_t &) const;
		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
path_step cath::test::return_path_matrix_test_suite_fixture::make_random_step() const {
	auto       rng           = default_random_engine{ random_device{}() };
	const auto random_triple = uniform_int_distribution<size_t>{ 0, 2 }( rng );
	switch ( random_triple ) {
		case ( 0 ) : { return path_step::ALIGN_PAIR;         }
		case ( 1 ) : { return path_step::INSERT_INTO_FIRST;  }
		default    : { return path_step::INSERT_INTO_SECOND; }
	}
}

/// \brief TODOCUMENT
return_path_matrix cath::test::return_path_matrix_test_suite_fixture::make_random_return_path_matrix(const size_t &arg_length_a,    ///< TODOCUMENT
                                                                                                     const size_t &arg_length_b,    ///< TODOCUMENT
                                                                                                     const size_t &arg_window_width ///< TODOCUMENT
                                                                                                     ) const {
	return_path_matrix new_path(arg_length_a, arg_length_b, arg_window_width);
	for (size_t ctr_a = 0; ctr_a < arg_length_a; ++ctr_a) {
		const size_size_pair b_window_start_and_stop = get_b_window_start_and_stop_for_a_index(
			new_path,
			ctr_a
		);
		const size_t &start = b_window_start_and_stop.first;
		const size_t &stop  = b_window_start_and_stop.second;
		for (size_t ctr_b = start; ctr_b <= stop; ++ctr_b) {
			path_step possible_random_step = make_random_step();
			while ( ( possible_random_step == path_step::INSERT_INTO_FIRST && ctr_b == start ) || ( possible_random_step == path_step::INSERT_INTO_SECOND && ctr_b == stop ) ) {
				possible_random_step = make_random_step();
			}
			new_path.set_path_step_towards_end_at_point(ctr_a, ctr_b, possible_random_step);
		}
	}
	return new_path;
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(return_path_matrix_test_suite, cath::test::return_path_matrix_test_suite_fixture)

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(basic) {
//	for (size_t ctr = 0; ctr < 5; ++ctr) {
//		cerr << make_random_return_path_matrix(13,  9, 7) << endl;
//		cerr << endl << endl << endl << endl << endl << endl;
//		cerr << make_random_return_path_matrix( 9, 13, 7) << endl;
//		cerr << endl << endl << endl << endl << endl << endl;
//	}
//}

/// \brief Important test that captured previous strange behaviour of attempting to implement
///        max_path_step_score with Boost Lambda (including bind()) in min_element
BOOST_AUTO_TEST_CASE(max_score_in_path_step_score_map_ones) {
	const path_step_score_map scores_of_path_step = { { path_step::ALIGN_PAIR,          1 },
	                                                  { path_step::INSERT_INTO_FIRST,  -1 },
	                                                  { path_step::INSERT_INTO_SECOND, -1 } };
	BOOST_CHECK_EQUAL( 1, max_path_step_score(scores_of_path_step) );
}

/// \brief Important test that captured previous strange behaviour of attempting to implement
///        max_path_step_score with Boost Lambda (including bind()) in min_element
BOOST_AUTO_TEST_CASE(max_score_in_path_step_score_map_real) {
	const path_step_score_map scores_of_path_step = { { path_step::ALIGN_PAIR,            0 },
	                                                  { path_step::INSERT_INTO_FIRST,  1564 },
	                                                  { path_step::INSERT_INTO_SECOND,  -50 } };
	BOOST_CHECK_EQUAL( 1564, max_path_step_score(scores_of_path_step) );
}

BOOST_AUTO_TEST_SUITE_END()
