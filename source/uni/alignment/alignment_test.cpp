/// \file
/// \brief The alignment test suite

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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

#include "alignment/io/align_scaffold.hpp"
#include "alignment/pair_alignment.hpp"
#include "alignment/test/alignment_fixture.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/pair_insertion_operator.hpp"
#include "common/size_t_literal.hpp"
#include "common/type_aliases.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "test/boost_addenda/boost_check_equal_ranges.hpp"
#include "test/boost_addenda/boost_check_no_throw_diag.hpp"
#include "test/test_tools.hpp"

#include <utility>

using namespace boost::test_tools;
using namespace cath;
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
			                                const size_size_pair_opt &);

			/// \brief Check that the has_position... and get_position methods work for an alignment by comparing to the lists from which it was made
			void check_details_match_lists(const alignment            &arg_aln,
			                               const aln_posn_opt_vec_vec &arg_lists
			                               ) {
				// Require that the sizes all match the length of the alignment
				for (const aln_posn_opt_vec &list : arg_lists) {
					BOOST_REQUIRE_EQUAL( arg_aln.length(), list.size() );
				}

				// Loop through the length of the alignment
				for (const size_t &index_ctr : indices( arg_aln.length() ) ) {
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
					for (const size_t &entry_ctr : indices( arg_lists.size() ) ) {
						// Check has_position_of_index_of_entry()
						const aln_posn_opt entry = arg_lists[entry_ctr][index_ctr];
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

	}  // namespace test
}  // namespace cath

/// \brief Check the consecutive position functions get the same result for the specified alignment as specified
void cath::test::alignment_test_suite_fixture::check_consecutive_position(const alignment          &arg_alignment, ///< The alignment to search
                                                                          const size_size_pair_opt &arg_expected   ///< The correct result
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
	BOOST_CHECK( true );
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
	for (const size_t &ctr : indices( aln_list_a.size() ) ) {
		const aln_posn_opt &a_posn = aln_list_a[ctr];
		const aln_posn_opt &b_posn = aln_list_b[ctr];

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
		for (const size_t &aln_ctr : indices( example_scores.size() ) ) {
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

BOOST_AUTO_TEST_SUITE(alignment_functions_test_suite)

BOOST_AUTO_TEST_CASE(has_positions_of_entry_in_index_range_works) {
	const auto the_aln = alignment_of_scaffold_lines( {
		"  XXXX",
		"XX  XX"
	} );

	BOOST_CHECK_THROW( has_positions_of_entry_in_index_range( the_aln, 0, 1, 0 ), invalid_argument_exception );
	BOOST_CHECK_THROW( has_positions_of_entry_in_index_range( the_aln, 1, 5, 4 ), invalid_argument_exception );

	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 0, 0, 0 ), false );
	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 1, 1, 1 ), false );

	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 0, 0, 2 ), false );
	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 0, 0, 3 ), true  );

	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 1, 1, 4 ), true  );
	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 1, 2, 4 ), false );
	BOOST_CHECK_EQUAL( has_positions_of_entry_in_index_range( the_aln, 1, 2, 5 ), true  );
}

BOOST_AUTO_TEST_CASE(entries_present_at_index_works) {
	const auto the_aln = alignment_of_scaffold_lines( {
		"  XXXX",
		"XX  XX"
	} );

	BOOST_CHECK_EQUAL_RANGES( entries_present_at_index( the_aln, 0 ), size_vec{    1 } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_at_index( the_aln, 1 ), size_vec{    1 } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_at_index( the_aln, 2 ), size_vec{ 0    } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_at_index( the_aln, 3 ), size_vec{ 0    } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_at_index( the_aln, 4 ), size_vec{ 0, 1 } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_at_index( the_aln, 5 ), size_vec{ 0, 1 } );
}

BOOST_AUTO_TEST_CASE(entries_present_in_index_range_works) {
	const auto the_aln = alignment_of_scaffold_lines( {
		"  XXXX",
		"XX  XX"
	} );

	BOOST_CHECK_THROW( entries_present_in_index_range( the_aln, 1, 0 ), invalid_argument_exception );
	BOOST_CHECK_THROW( entries_present_in_index_range( the_aln, 5, 4 ), invalid_argument_exception );

	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 0, 0 ), size_vec{      } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 0, 2 ), size_vec{    1 } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 0, 3 ), size_vec{ 0, 1 } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 1, 1 ), size_vec{      } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 1, 4 ), size_vec{ 0, 1 } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 2, 4 ), size_vec{ 0    } );
	BOOST_CHECK_EQUAL_RANGES( entries_present_in_index_range( the_aln, 2, 5 ), size_vec{ 0, 1 } );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
