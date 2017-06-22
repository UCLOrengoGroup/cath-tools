/// \file
/// \brief The seq_seg_run test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/algorithm/is_uniq_for_unordered.hpp"
#include "seq/seq_seg_run.hpp"

namespace cath { namespace test { } }

using namespace cath::common;
using namespace cath::seq;
using namespace cath::test;

using boost::lexical_cast;
using std::pair;
using std::string;

namespace cath {
	namespace seq {
		using seq_seg_seq_seg_pair = pair< seq_seg, seq_seg >;
	}
}

namespace cath {
	namespace test {

		/// \brief The seq_seg_run_test_suite_fixture to assist in testing seq_seg_run
		struct seq_seg_run_test_suite_fixture {
		protected:
			~seq_seg_run_test_suite_fixture() noexcept = default;

			static constexpr seq_seg seg_zero_one   { 0, 1 };
			static constexpr seq_seg seg_zero_two   { 0, 2 };
			static constexpr seq_seg seg_zero_three { 0, 3 };
			static constexpr seq_seg seg_one_two    { 1, 2 };
			static constexpr seq_seg seg_two_three  { 2, 3 };
			static constexpr seq_seg seg_three_four { 3, 4 };

			/// \brief Example of two single segments: sngls_separated_1_before
			// |--   |
			// |   --|
			static constexpr seq_seg_seq_seg_pair sngls_separated_1_before     { seg_zero_one,   seg_three_four };

			/// \brief Example of two single segments: sngls_separated_2_before
			// |   --|
			// |--   |
			static constexpr seq_seg_seq_seg_pair sngls_separated_2_before     { seg_three_four, seg_zero_one   };

			/// \brief Example of two single segments: sngls_touching_1_before
			// |--   |
			// |  -- |
			static constexpr seq_seg_seq_seg_pair sngls_touching_1_before      { seg_zero_one,   seg_two_three  };

			/// \brief Example of two single segments: sngls_touching_2_before
			// |  -- |
			// |--   |
			static constexpr seq_seg_seq_seg_pair sngls_touching_2_before      { seg_two_three,  seg_zero_one   };

			/// \brief Example of two single segments: sngls_ovrlp_1_before
			// |--   |
			// | --  |
			static constexpr seq_seg_seq_seg_pair sngls_ovrlp_1_before         { seg_zero_one,   seg_one_two    };

			/// \brief Example of two single segments: sngls_ovrlp_2_before
			// | --  |
			// |--   |
			static constexpr seq_seg_seq_seg_pair sngls_ovrlp_2_before         { seg_one_two,    seg_zero_one   };

			/// \brief Example of two single segments: sngls_ovrhng_1_starts_before
			// |---  |
			// | --  |
			static constexpr seq_seg_seq_seg_pair sngls_ovrhng_1_starts_before { seg_zero_two,   seg_one_two    };

			/// \brief Example of two single segments: sngls_ovrhng_2_starts_before
			// | --  |
			// |---  |
			static constexpr seq_seg_seq_seg_pair sngls_ovrhng_2_starts_before { seg_one_two,    seg_zero_two   };

			/// \brief Example of two single segments: sngls_ovrhng_2_stops_after
			// |--   |
			// |---  |
			static constexpr seq_seg_seq_seg_pair sngls_ovrhng_2_stops_after   { seg_zero_one,   seg_zero_two   };

			/// \brief Example of two single segments: sngls_ovrhng_1_stops_after
			// |---  |
			// |--   |
			static constexpr seq_seg_seq_seg_pair sngls_ovrhng_1_stops_after   { seg_zero_two,   seg_zero_one   };

			/// \brief Example of two single segments: sngls_nested_1_within
			// |---- |
			// | --  |
			static constexpr seq_seg_seq_seg_pair sngls_nested_1_within        { seg_zero_three, seg_one_two    };

			/// \brief Example of two single segments: sngls_nested_2_within
			// | --  |
			// |---- |
			static constexpr seq_seg_seq_seg_pair sngls_nested_2_within        { seg_one_two,    seg_zero_three };

			/// \brief Example of two single segments: sngls_same
			// |--   |
			// |--   |
			static constexpr seq_seg_seq_seg_pair sngls_same                   { seg_zero_one,   seg_zero_one   };

			/// \brief Make a seq_seg_run of the first seq_seg in the specified pair
			seq_seg_run sngl_run_1(const seq_seg_seq_seg_pair &arg_seq_segs ///< The pair of single seq_segs to query
			                       ) {
				return { { arg_seq_segs.first } };
			}

			/// \brief Make a seq_seg_run of the second seq_seg in the specified pair
			seq_seg_run sngl_run_2(const seq_seg_seq_seg_pair &arg_seq_segs ///< The pair of single seq_segs to query
			                      ) {
				return { { arg_seq_segs.second } };
			}
		};

	} // namespace test
} // namespace cath

constexpr seq_seg              seq_seg_run_test_suite_fixture::seg_zero_one;
constexpr seq_seg              seq_seg_run_test_suite_fixture::seg_zero_two;
constexpr seq_seg              seq_seg_run_test_suite_fixture::seg_zero_three;
constexpr seq_seg              seq_seg_run_test_suite_fixture::seg_one_two;
constexpr seq_seg              seq_seg_run_test_suite_fixture::seg_two_three;
constexpr seq_seg              seq_seg_run_test_suite_fixture::seg_three_four;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_separated_1_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_separated_2_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_touching_1_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_touching_2_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_ovrlp_1_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_ovrlp_2_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_ovrhng_1_starts_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_ovrhng_2_starts_before;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_ovrhng_2_stops_after;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_ovrhng_1_stops_after;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_nested_1_within;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_nested_2_within;
constexpr seq_seg_seq_seg_pair seq_seg_run_test_suite_fixture::sngls_same;

BOOST_FIXTURE_TEST_SUITE(seq_seg_run_test_suite, seq_seg_run_test_suite_fixture)


BOOST_AUTO_TEST_CASE(make_seq_seg_run_from_res_indices_works) {
	const seq_seg_run a{ seq_seg_vec{ {  0, 1 }, { 4, 6 } } };
	const seq_seg_run b{ seq_seg_vec{ {  3, 5 }           } };
	BOOST_CHECK_EQUAL( make_seq_seg_run_from_res_indices( 0, 1, 4, 6 ), a );
	BOOST_CHECK_EQUAL( make_seq_seg_run_from_res_indices( 3,       5 ), b );
}


BOOST_AUTO_TEST_CASE(converts_to_string_correctly) {
	const seq_seg_run a{ seq_seg_vec{ {  100,  199 }, {  300,  399 } } };
	BOOST_CHECK_EQUAL( to_string           ( a ), "seq_seg_run[100-199,300-399]" );
	BOOST_CHECK_EQUAL( lexical_cast<string>( a ), "seq_seg_run[100-199,300-399]" );
}

BOOST_AUTO_TEST_CASE(basic) {
	const seq_seg_run a{ seq_seg_vec{ {  100,  199 }, {  300,  399 } } };
	const seq_seg_run b{ seq_seg_vec{ {  190,  209 }, {  390,  409 } } };
	const seq_seg_run c{ seq_seg_vec{ {  200,  299 }, {  400,  499 } } };
	const seq_seg_run d{ seq_seg_vec{ { 1100, 1199 }, { 1400, 1499 } } };

	BOOST_CHECK      (   are_overlapping               ( a, a )      );
	BOOST_CHECK      (   are_overlapping               ( a, b )      );
	BOOST_CHECK      ( ! are_overlapping               ( a, c )      );
	BOOST_CHECK      ( ! are_overlapping               ( a, d )      );
	BOOST_CHECK      (   are_overlapping               ( b, a )      );
	BOOST_CHECK      (   are_overlapping               ( b, b )      );
	BOOST_CHECK      (   are_overlapping               ( b, c )      );
	BOOST_CHECK      ( ! are_overlapping               ( b, d )      );
	BOOST_CHECK      ( ! are_overlapping               ( c, a )      );
	BOOST_CHECK      (   are_overlapping               ( c, b )      );
	BOOST_CHECK      (   are_overlapping               ( c, c )      );
	BOOST_CHECK      ( ! are_overlapping               ( c, d )      );
	BOOST_CHECK      ( ! are_overlapping               ( d, a )      );
	BOOST_CHECK      ( ! are_overlapping               ( d, b )      );
	BOOST_CHECK      ( ! are_overlapping               ( d, c )      );
	BOOST_CHECK      (   are_overlapping               ( d, d )      );

	BOOST_CHECK_EQUAL(   overlap_by                    ( a, a ), 200 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( a, b ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( a, c ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( a, d ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, a ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, b ),  40 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, c ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( b, d ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, a ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, b ),  20 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, c ), 200 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( c, d ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, a ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, b ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, c ),   0 );
	BOOST_CHECK_EQUAL(   overlap_by                    ( d, d ), 200 );

	BOOST_CHECK_EQUAL(   shorter_length                ( a, a ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( a, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( a, c ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( a, d ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, a ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, c ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( b, d ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, a ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, c ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( c, d ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, a ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, b ),  40 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, c ), 200 );
	BOOST_CHECK_EQUAL(   shorter_length                ( d, d ), 200 );

	BOOST_CHECK_EQUAL(   longer_length                 ( a, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( a, b ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( a, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( a, d ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, b ),  40 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( b, d ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, b ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( c, d ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, a ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, b ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, c ), 200 );
	BOOST_CHECK_EQUAL(   longer_length                 ( d, d ), 200 );

	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, a ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, b ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( a, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, a ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, b ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, c ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( b, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, b ), 0.5 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, c ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( c, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, b ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_shorter ( d, d ), 1.0 );

	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, a ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, b ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( a, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, a ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, b ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, c ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( b, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, b ), 0.1 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, c ), 1.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( c, d ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, a ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, b ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, c ), 0.0 );
	BOOST_CHECK_EQUAL(   fraction_overlap_over_longer  ( d, d ), 1.0 );
}

BOOST_AUTO_TEST_CASE(problem_case) {
	BOOST_CHECK(
		are_overlapping(
			make_seq_seg_run_from_res_indices( 0, 1, 4, 6 ),
			make_seq_seg_run_from_res_indices( 3,       5 )
		)
	);
}

BOOST_AUTO_TEST_SUITE(covering_functions)


BOOST_AUTO_TEST_SUITE(singles)

// |--   |
// |   --|
BOOST_AUTO_TEST_CASE(work_on_singles_separated_1_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_separated_1_before     ), sngl_run_2( sngls_separated_1_before     ) ) );
	BOOST_CHECK( ! one_covers_other                  ( sngl_run_1( sngls_separated_1_before     ), sngl_run_2( sngls_separated_1_before     ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_separated_1_before     ), sngl_run_2( sngls_separated_1_before     ) ) );
}

// |   --|
// |--   |
BOOST_AUTO_TEST_CASE(work_on_singles_separated_2_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_separated_2_before     ), sngl_run_2( sngls_separated_2_before     ) ) );
	BOOST_CHECK( ! one_covers_other                  ( sngl_run_1( sngls_separated_2_before     ), sngl_run_2( sngls_separated_2_before     ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_separated_2_before     ), sngl_run_2( sngls_separated_2_before     ) ) );
}

// |--   |
// |  -- |
BOOST_AUTO_TEST_CASE(work_on_singles_touching_1_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_touching_1_before      ), sngl_run_2( sngls_touching_1_before      ) ) );
	BOOST_CHECK( ! one_covers_other                  ( sngl_run_1( sngls_touching_1_before      ), sngl_run_2( sngls_touching_1_before      ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_touching_1_before      ), sngl_run_2( sngls_touching_1_before      ) ) );
}

// |  -- |
// |--   |
BOOST_AUTO_TEST_CASE(work_on_singles_touching_2_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_touching_2_before      ), sngl_run_2( sngls_touching_2_before      ) ) );
	BOOST_CHECK( ! one_covers_other                  ( sngl_run_1( sngls_touching_2_before      ), sngl_run_2( sngls_touching_2_before      ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_touching_2_before      ), sngl_run_2( sngls_touching_2_before      ) ) );
}

// |--   |
// | --  |
BOOST_AUTO_TEST_CASE(work_on_singles_overlapping_1_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_ovrlp_1_before         ), sngl_run_2( sngls_ovrlp_1_before         ) ) );
	BOOST_CHECK( ! one_covers_other                  ( sngl_run_1( sngls_ovrlp_1_before         ), sngl_run_2( sngls_ovrlp_1_before         ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_ovrlp_1_before         ), sngl_run_2( sngls_ovrlp_1_before         ) ) );
}

// | --  |
// |--   |
BOOST_AUTO_TEST_CASE(work_on_singles_overlapping_2_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_ovrlp_2_before         ), sngl_run_2( sngls_ovrlp_2_before         ) ) );
	BOOST_CHECK( ! one_covers_other                  ( sngl_run_1( sngls_ovrlp_2_before         ), sngl_run_2( sngls_ovrlp_2_before         ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_ovrlp_2_before         ), sngl_run_2( sngls_ovrlp_2_before         ) ) );
}

// |---  |
// | --  |
BOOST_AUTO_TEST_CASE(work_on_singles_overhang_1_starts_before) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_ovrhng_1_starts_before ), sngl_run_2( sngls_ovrhng_1_starts_before ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_ovrhng_1_starts_before ), sngl_run_2( sngls_ovrhng_1_starts_before ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_ovrhng_1_starts_before ), sngl_run_2( sngls_ovrhng_1_starts_before ) ) );
}

// | --  |
// |---  |
BOOST_AUTO_TEST_CASE(work_on_singles_overhang_2_starts_before) {
	BOOST_CHECK(   first_is_not_outside_second       ( sngl_run_1( sngls_ovrhng_2_starts_before ), sngl_run_2( sngls_ovrhng_2_starts_before ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_ovrhng_2_starts_before ), sngl_run_2( sngls_ovrhng_2_starts_before ) ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( sngl_run_1( sngls_ovrhng_2_starts_before ), sngl_run_2( sngls_ovrhng_2_starts_before ) ) );
}

// |--   |
// |---  |
BOOST_AUTO_TEST_CASE(work_on_singles_overhang_2_stops_after) {
	BOOST_CHECK(   first_is_not_outside_second       ( sngl_run_1( sngls_ovrhng_2_stops_after   ), sngl_run_2( sngls_ovrhng_2_stops_after   ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_ovrhng_2_stops_after   ), sngl_run_2( sngls_ovrhng_2_stops_after   ) ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( sngl_run_1( sngls_ovrhng_2_stops_after   ), sngl_run_2( sngls_ovrhng_2_stops_after   ) ) );
}

// |---  |
// |--   |
BOOST_AUTO_TEST_CASE(work_on_singles_overhang_1_stops_after) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_ovrhng_1_stops_after   ), sngl_run_2( sngls_ovrhng_1_stops_after   ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_ovrhng_1_stops_after   ), sngl_run_2( sngls_ovrhng_1_stops_after   ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_ovrhng_1_stops_after   ), sngl_run_2( sngls_ovrhng_1_stops_after   ) ) );
}

// |---- |
// | --  |
BOOST_AUTO_TEST_CASE(work_on_singles_nested_1_within) {
	BOOST_CHECK( ! first_is_not_outside_second       ( sngl_run_1( sngls_nested_1_within        ), sngl_run_2( sngls_nested_1_within        ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_nested_1_within        ), sngl_run_2( sngls_nested_1_within        ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_nested_1_within        ), sngl_run_2( sngls_nested_1_within        ) ) );
}

// | --  |
// |---- |
BOOST_AUTO_TEST_CASE(work_on_singles_nested_2_within) {
	BOOST_CHECK(   first_is_not_outside_second       ( sngl_run_1( sngls_nested_2_within        ), sngl_run_2( sngls_nested_2_within        ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_nested_2_within        ), sngl_run_2( sngls_nested_2_within        ) ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( sngl_run_1( sngls_nested_2_within        ), sngl_run_2( sngls_nested_2_within        ) ) );
}

// |--   |
// |--   |
BOOST_AUTO_TEST_CASE(work_on_singles_identical) {
	BOOST_CHECK(   first_is_not_outside_second       ( sngl_run_1( sngls_same              ), sngl_run_2( sngls_same              ) ) );
	BOOST_CHECK(   one_covers_other                  ( sngl_run_1( sngls_same              ), sngl_run_2( sngls_same              ) ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( sngl_run_1( sngls_same              ), sngl_run_2( sngls_same              ) ) );
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(multis)


BOOST_AUTO_TEST_CASE(same) {
	const seq_seg_run a{ seq_seg_vec{ {  100,  199 }, {  300,  399 } } };
	BOOST_CHECK(   first_is_not_outside_second       ( a, a ) );
	BOOST_CHECK(   one_covers_other                  ( a, a ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( a, a ) );
}

BOOST_AUTO_TEST_CASE(extra_seg_at_start) {
	const seq_seg_run a{ seq_seg_vec{ { 30, 59 }, {  100,  199 }, {  300,  399 }               } };
	const seq_seg_run b{ seq_seg_vec{             {  100,  199 }, {  300,  399 }               } };
	BOOST_CHECK( ! first_is_not_outside_second       ( a, b ) );
	BOOST_CHECK(   one_covers_other                  ( a, b ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( a, b ) );

	BOOST_CHECK(   first_is_not_outside_second       ( b, a ) );
	BOOST_CHECK(   one_covers_other                  ( b, a ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( b, a ) );
}

BOOST_AUTO_TEST_CASE(extra_seg_at_end) {
	const seq_seg_run a{ seq_seg_vec{             {  100,  199 }, {  300,  399 }, { 500, 599 } } };
	const seq_seg_run b{ seq_seg_vec{             {  100,  199 }, {  300,  399 }               } };
	BOOST_CHECK( ! first_is_not_outside_second       ( a, b ) );
	BOOST_CHECK(   one_covers_other                  ( a, b ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( a, b ) );

	BOOST_CHECK(   first_is_not_outside_second       ( b, a ) );
	BOOST_CHECK(   one_covers_other                  ( b, a ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( b, a ) );
}



BOOST_AUTO_TEST_CASE(gap_start_offset) {
	const seq_seg_run a{ seq_seg_vec{ {  100,  198 }, {  300,  399 } } };
	const seq_seg_run b{ seq_seg_vec{ {  100,  199 }, {  300,  399 } } };

	BOOST_CHECK(   first_is_not_outside_second       ( a, b ) );
	BOOST_CHECK(   one_covers_other                  ( a, b ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( a, b ) );

	BOOST_CHECK( ! first_is_not_outside_second       ( b, a ) );
	BOOST_CHECK(   one_covers_other                  ( b, a ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( b, a ) );
}

BOOST_AUTO_TEST_CASE(gap_stop_offset) {
	const seq_seg_run a{ seq_seg_vec{ {  100,  199 }, {  301,  399 } } };
	const seq_seg_run b{ seq_seg_vec{ {  100,  199 }, {  300,  399 } } };

	BOOST_CHECK(   first_is_not_outside_second       ( a, b ) );
	BOOST_CHECK(   one_covers_other                  ( a, b ) );
	BOOST_CHECK(   first_is_shorter_and_within_second( a, b ) );

	BOOST_CHECK( ! first_is_not_outside_second       ( b, a ) );
	BOOST_CHECK(   one_covers_other                  ( b, a ) );
	BOOST_CHECK( ! first_is_shorter_and_within_second( b, a ) );
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(multis)

BOOST_AUTO_TEST_CASE(basic) {
	const auto hashes = { calc_hash( make_seq_seg_run_from_res_indices( 100,             399 ) ),
	                      calc_hash( make_seq_seg_run_from_res_indices( 101,             399 ) ),
	                      calc_hash( make_seq_seg_run_from_res_indices( 100,             398 ) ),
	                      calc_hash( make_seq_seg_run_from_res_indices( 100,  199, 300,  399 ) ),
	                      calc_hash( make_seq_seg_run_from_res_indices( 100,  198, 300,  399 ) ),
	                      calc_hash( make_seq_seg_run_from_res_indices( 100,  199, 301,  399 ) ) };

	BOOST_CHECK( is_uniq_for_unordered( hashes ) );
}

BOOST_AUTO_TEST_SUITE_END()



// BOOST_AUTO_TEST_SUITE(middle)

BOOST_AUTO_TEST_CASE(middle) {
	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 2, 2       )         ),  2.0 );
	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 2, 3       )         ),  2.5 );
	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 1, 3       )         ),  2.0 );

	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 1,       9 )         ),  5.0 );
	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 1, 3, 7, 9 )         ),  5.0 );

	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 1,       9, 21, 29 ) ), 15.0 );
	BOOST_CHECK_EQUAL( middle_index( make_seq_seg_run_from_res_indices( 1, 3, 7, 9, 21, 29 ) ), 17.0 );
}

// BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE_END()
