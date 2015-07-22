/// \file
/// \brief The alignment test suite

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/test/output_test_stream.hpp>

#include "alignment/pair_alignment.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/pair_insertion_operator.h"
#include "common/size_t_literal.h"
#include "common/test_tools.h"
#include "exception/invalid_argument_exception.h"
#include "test/alignment_fixture.h"

#include <utility>

using namespace boost::test_tools;
using namespace cath::align;
using namespace cath::common;
using namespace cath::common::test;
using namespace std;

using boost::none;
using boost::numeric_cast;

namespace cath {
	namespace test {

		/// \brief The alignment_test_suite_fixture to assist in testing alignment
		struct alignment_test_suite_fixture : protected alignment_fixture {
			~alignment_test_suite_fixture() noexcept = default;

			void check_consecutive_position(const alignment &,
			                                const opt_size_size_pair &);

			/// \brief Check that the has_position... and get_position methods work for an alignment by comparing to the lists from which it was made
			void check_details_match_lists(const alignment            &arg_aln,
			                               const opt_aln_posn_vec_vec &arg_lists
			                               ) {
				// Require that the sizes all match the length of the alignment
				for (const opt_aln_posn_vec &list : arg_lists) {
					BOOST_REQUIRE_EQUAL( arg_aln.length(), list.size() );
				}

				// Loop through the length of the alignment
				for (size_t index_ctr = 0; index_ctr < arg_aln.length(); ++index_ctr) {
					// If this is a pair alignment, check has_both_positions_of_index()
					const bool is_last  = (index_ctr +1 == arg_aln.length());
					const bool is_pair_alignment = (arg_lists.size() == 2_z);
					if ( is_pair_alignment ) {
						BOOST_CHECK_EQUAL(
							( arg_lists[alignment::PAIR_A_IDX][index_ctr] && arg_lists[alignment::PAIR_B_IDX][index_ctr] ),
							has_both_positions_of_index(arg_aln, index_ctr )
						);
					}

					// Loop over the entries
					for (size_t entry_ctr = 0; entry_ctr < arg_lists.size(); ++entry_ctr) {
						// Check has_position_of_index_of_entry()
						const opt_aln_posn entry = arg_lists[entry_ctr][index_ctr];
						BOOST_REQUIRE_EQUAL(
							static_cast<bool>( entry ),
							has_position_of_entry_of_index( arg_aln, entry_ctr, index_ctr )
						);
						if ( is_last ) {
							BOOST_REQUIRE_EQUAL(
								static_cast<bool>( entry ),
								has_last_position_of_entry(arg_aln, entry_ctr)
							);
						}

						// If has_entry, check the get_position_of_index_of_entry() result
						if ( entry ) {
							BOOST_CHECK_EQUAL(
								*entry,
								get_position_of_entry_of_index( arg_aln, entry_ctr, index_ctr )
							);
							if (is_last) {
								BOOST_CHECK_EQUAL(
									*entry,
									get_last_position_of_entry(arg_aln, entry_ctr)
								);
							}
						}
						// Otherwise check get_position_of_index_of_entry() throws
						else {
		//					BOOST_CHECK_THROW(arg_aln.get_position_of_entry_of_index(entry_ctr, index_ctr), invalid_argument_exception);
							BOOST_CHECK_THROW( get_position_of_entry_of_index( arg_aln, entry_ctr, index_ctr ), invalid_argument_exception );
							if ( is_last ) {
								BOOST_CHECK_THROW(
									get_last_position_of_entry(arg_aln, entry_ctr),
									invalid_argument_exception
								);
							}
						}

						// If this is a pair alignment, check has_a_position_of_index() and get_a_position_of_index()
						if ( is_pair_alignment ) {
							BOOST_REQUIRE_EQUAL(
								static_cast<bool>( entry ),
								(entry_ctr == alignment::PAIR_A_IDX) ? has_a_position_of_index(arg_aln, index_ctr) : has_b_position_of_index(arg_aln, index_ctr)
							);
							if (is_last) {
								BOOST_REQUIRE_EQUAL(
									static_cast<bool>( entry ),
									(entry_ctr == alignment::PAIR_A_IDX) ? has_last_a_position(arg_aln) : has_last_b_position(arg_aln)
								);
							}

							if ( entry ) {
								BOOST_CHECK_EQUAL(
									*entry,
									(entry_ctr == alignment::PAIR_A_IDX) ? get_a_position_of_index(arg_aln, index_ctr) : get_b_position_of_index(arg_aln, index_ctr)
								);
								if (is_last) {
									BOOST_REQUIRE_EQUAL(
										*entry,
										(entry_ctr == alignment::PAIR_A_IDX) ? get_last_a_position(arg_aln) : get_last_b_position(arg_aln)
									);
								}
							}
							else {
								BOOST_CHECK_THROW(
									(entry_ctr == alignment::PAIR_A_IDX) ? get_a_position_of_index( arg_aln, index_ctr ) : get_b_position_of_index( arg_aln, index_ctr ),
									invalid_argument_exception
								);
							}
						}
					}
				}
			}
		};

	}
}

/// \brief Check the consecutive position functions get the same result for the specified alignment as specified
void cath::test::alignment_test_suite_fixture::check_consecutive_position(const alignment          &arg_alignment, ///< The alignment to search
                                                                          const opt_size_size_pair &arg_expected   ///< The correct result
                                                                          ) {

	if ( arg_expected ) {
		BOOST_CHECK_EQUAL( first_non_consecutive_entry_positions( arg_alignment ), arg_expected               );
		BOOST_CHECK_THROW( check_entry_positions_are_consecutive( arg_alignment ), invalid_argument_exception );
	}
	else {
		BOOST_CHECK_EQUAL        ( first_non_consecutive_entry_positions( arg_alignment ), none );
		BOOST_CHECK_NO_THROW_DIAG( check_entry_positions_are_consecutive( arg_alignment )       );
	}
}

BOOST_FIXTURE_TEST_SUITE(alignment_test_suite, cath::test::alignment_test_suite_fixture)

/// \brief Just test a trivially true equality to ensure that all the fixture's setup works OK
BOOST_AUTO_TEST_CASE(fixture_constructions_work) {
	BOOST_CHECK_EQUAL(1, 1);
}

/// \brief Check that the list constructor throws if passed lists of differing lengths
BOOST_AUTO_TEST_CASE(different_size_throws) {
	BOOST_CHECK_THROW( alignment test_aln( { aln_list_a  , aln_list_long } ), invalid_argument_exception );
	BOOST_CHECK_THROW( alignment test_aln( { aln_list_b  , aln_list_long } ), invalid_argument_exception );
	BOOST_CHECK_THROW( alignment test_aln( {aln_list_long, aln_list_a    } ), invalid_argument_exception );
	BOOST_CHECK_THROW( alignment test_aln( {aln_list_long, aln_list_b    } ), invalid_argument_exception );
}

/// \brief Check that the size method works as expected for pair alignments
BOOST_AUTO_TEST_CASE(size) {
	BOOST_CHECK_EQUAL( aln_list_a.size(), aln_a_a.length() );

	BOOST_CHECK_EQUAL( aln_list_a.size(), aln_a_b.length() );
	BOOST_CHECK_EQUAL( aln_list_b.size(), aln_a_b.length() );

	BOOST_CHECK_EQUAL( aln_list_a.size(), aln_b_a.length() );
	BOOST_CHECK_EQUAL( aln_list_b.size(), aln_b_a.length() );
}

/// \brief Check that the equality operator works as expected for pair alignments
BOOST_AUTO_TEST_CASE(equality) {
	const auto alignments = { aln_a_a, aln_a_b, aln_b_a, aln_long_long };
	check_equality_operators_on_diff_vals_range( alignments );
}

/// \brief Check that the three append_operators behave as expected for a pair alignment
BOOST_AUTO_TEST_CASE(append_a_b_and_both) {
	BOOST_REQUIRE_EQUAL(aln_list_a.size(), aln_list_b.size());
	alignment new_aln_a_b(alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT);
	alignment new_aln_b_a(alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT);
	for ( size_t ctr = 0; ctr < aln_list_a.size(); ++ctr ) {
		const opt_aln_posn &a_posn = aln_list_a[ctr];
		const opt_aln_posn &b_posn = aln_list_b[ctr];

		if ( a_posn  && b_posn ) {
			append_position_both( new_aln_a_b, *a_posn, *b_posn );
			append_position_both( new_aln_b_a, *b_posn, *a_posn );
		}
		else if ( a_posn ) {
			append_position_a( new_aln_a_b, *a_posn );
			append_position_b( new_aln_b_a, *a_posn );
		}
		else if ( b_posn ) {
			append_position_b( new_aln_a_b, *b_posn );
			append_position_a( new_aln_b_a, *b_posn );
		}
		else {
			BOOST_REQUIRE(false); // If we hit this line, something has gone badly wrong
		}
	}

	BOOST_CHECK_EQUAL( aln_a_b, new_aln_a_b );
	BOOST_CHECK_EQUAL( aln_b_a, new_aln_b_a );
}

/// \brief Check that the three append_operators for a pair alignment
BOOST_AUTO_TEST_CASE(has_position_and_get_position_a_a) {
	check_details_match_lists(aln_a_a, { aln_list_a, aln_list_a } );
}
BOOST_AUTO_TEST_CASE(has_position_and_get_position_a_b) {
	check_details_match_lists(aln_a_b, { aln_list_a, aln_list_b } );
}
BOOST_AUTO_TEST_CASE(has_position_and_get_position_b_a) {
	check_details_match_lists(aln_b_a, { aln_list_b, aln_list_a } );
}

/// \brief Check that the scores behave as expected
BOOST_AUTO_TEST_CASE(scores) {
	// Can't do aln_a_a because that requires 4 scored positions
	const vector<alignment> test_alignments = { aln_a_b, aln_b_a };
	for (const alignment &test_alignment : test_alignments) {
		BOOST_CHECK(!test_alignment.is_scored());
		alignment copy_alignment(test_alignment);
		BOOST_CHECK(!copy_alignment.is_scored());
		set_pair_alignment_duplicate_scores( copy_alignment, example_scores );
		BOOST_CHECK(copy_alignment.is_scored());
		size_t score_ctr = 0;
		for (size_t aln_ctr = 0; aln_ctr < example_scores.size(); ++aln_ctr) {
			if (has_both_positions_of_index(copy_alignment, aln_ctr)) {
				BOOST_CHECK_EQUAL(
					example_scores[score_ctr],
					get_score_of_entry_and_index( copy_alignment, alignment::PAIR_A_IDX, aln_ctr)
				);
				BOOST_CHECK_EQUAL(
					example_scores[score_ctr],
					get_score_of_entry_and_index( copy_alignment, alignment::PAIR_B_IDX, aln_ctr)
				);
				++score_ctr;
			}
		}
	}

	// This does not currently test set_non_raw_scores()
}

/// \brief Check that the insertion operator summarises alignments as expected
BOOST_AUTO_TEST_CASE(insertion_operator) {
	output_test_stream output;
	output << aln_a_a;
	BOOST_CHECK( output.is_equal( "alignment[4 positions: 0 <-> 0; 1 <-> 1; 2 <-> 2; 3 <-> 3]" ) );
	output << aln_a_b;
	BOOST_CHECK( output.is_equal( "alignment[4 positions: 0 <-> 0; 1 <-> 1; 2 <-> 2; 3 <-> ]"  ) );
	output << aln_b_a;
	BOOST_CHECK( output.is_equal( "alignment[4 positions: 0 <-> 0; 1 <-> 1; 2 <-> 2;  <-> 3]"  ) );
}

/// \brief Check that get_last_present_position_of_entry() works as expected
BOOST_AUTO_TEST_CASE( get_last_present_position_of_entry_works ) {
	BOOST_CHECK_EQUAL( 3_z, get_max_last_present_position( aln_a_a       ).get() );
	BOOST_CHECK_EQUAL( 3_z, get_max_last_present_position( aln_a_b       ).get() );
	BOOST_CHECK_EQUAL( 3_z, get_max_last_present_position( aln_b_a       ).get() );
	BOOST_CHECK_EQUAL( 5_z, get_max_last_present_position( aln_long_long ).get() );
}

/// \brief Check that first_non_consecutive_entry_positions works as expected
BOOST_AUTO_TEST_CASE( first_non_consecutive_entry_positions ) {
	check_consecutive_position(
		alignment( {
			{ 0_z, 1,    none },
			{ 0_z, none, 1    }
		} ),
		none
	);
	check_consecutive_position(
		alignment( {
			{ 0_z, 1 },
			{ 0_z, 2 }
		} ),
		make_pair( 1_z, 1_z )
	);
	check_consecutive_position(
		alignment( {
			{ 0_z, 2 },
			{ 0_z, 1 }
		} ),
		make_pair( 0_z, 1_z )
	);
	check_consecutive_position(
		alignment( {
			{ 5, 6 },
			{ 5, 6 }
		} ),
		make_pair( 0_z, 0_z )
	);
}

BOOST_AUTO_TEST_SUITE_END()