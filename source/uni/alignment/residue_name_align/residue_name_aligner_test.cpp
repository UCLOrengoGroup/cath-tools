/// \file


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

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/permutation.hpp>

#include "alignment/alignment.hpp"
#include "alignment/residue_name_align/residue_name_aligner.hpp"
#include "biocore/residue_name.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/size_t_literal.hpp"

#include <deque>
#include <vector>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::range::next_permutation;

namespace cath {
	namespace test {

		/// \brief The residue_name_aligner_test_suite_fixture to assist in testing residue_name_aligner
		struct residue_name_aligner_test_suite_fixture {
		private:
			/// \brief Checks that residue_name_aligner::residue_name_align() does the correct thing for all permutations of the residue lists.
			///
			/// This should be accessed via check_residue_name_aligner_results() or check_residue_name_aligner_throws()
			///
			/// This subroutine is mostly ready to handle more than two lists.
			///
			/// \todo modify the actual call to residue_name_aligner::residue_name_align() to handle multiple lists
			///       (after that subroutine has been altered to handle multiple lists)
			void do_check_residue_name_aligner(const residue_name_vec_vec  &arg_residue_lists,
			                                   const bool_deq_vec          &arg_correct_presence_lists,
			                                   const size_vec_vec          &arg_correct_answer_lists,
			                                   const bool                  &arg_should_throw
			                                   ) {
				const size_t num_lists = arg_residue_lists.size();
				BOOST_REQUIRE_EQUAL(num_lists, 2_z); /// This code isn't yet able to process more than two at a time
				if (!arg_should_throw) {
					BOOST_REQUIRE_EQUAL(num_lists, arg_correct_presence_lists.size());
					BOOST_REQUIRE_EQUAL(num_lists, arg_correct_answer_lists.size());
				}

				// Construct a vector containing the indices of the lists
				size_vec permutation_indices(num_lists, 0);
				for (const size_t &index_ctr : indices( num_lists ) ) {
					permutation_indices[index_ctr] = index_ctr;
				}

				// Loop over the permutations of the indices
				do {
					// If these residue lists should cause residue_name_aligner::residue_name_align() to throw then check they do
					if (arg_should_throw) {
						BOOST_CHECK_THROW(
							residue_name_aligner::residue_name_align(
								{ arg_residue_lists[ permutation_indices[ 0 ] ],
								  arg_residue_lists[ permutation_indices[ 1 ] ] }
							),
							invalid_argument_exception
						);
					}
					// Otherwise check the results from residue_name_aligner::residue_name_align()
					else {
						// Construct an alignment from this permutation of residue lists
						const alignment my_alignment = residue_name_aligner::residue_name_align(
							{ arg_residue_lists[permutation_indices[ 0 ]],
							  arg_residue_lists[permutation_indices[ 1 ]] }
						);
						const alignment::size_type num_positions = my_alignment.length();

						// Check each of the alignment entries in turn
						for (const size_t &index_ctr : indices( num_lists ) ) {
							// Grab the correct answer list under the current permutation
							const size_t permutation_index           = permutation_indices[index_ctr];
							const bool_deq &correct_presence_list = arg_correct_presence_lists[permutation_index];
							const size_vec &correct_answer_list      = arg_correct_answer_lists[permutation_index];

							const size_t correct_answer_size         = correct_answer_list.size();

							// Check that the number of positions match
							BOOST_CHECK_EQUAL(correct_answer_size, num_positions);

							// Check that each of the positions in the alignment match what is expected
							for (const size_t &position_ctr : indices( min( correct_answer_size, num_positions ) ) ) {
								const aln_posn_opt position     = my_alignment.position_of_entry_of_index( index_ctr, position_ctr );
								const bool         has_position = static_cast<bool>( position );
								BOOST_CHECK_EQUAL( correct_presence_list[position_ctr], has_position );

								if ( position ) {
									BOOST_CHECK_EQUAL(correct_answer_list[position_ctr], *position );
								}
							}
						}
					}
				} while ( next_permutation( permutation_indices ) );
			}

		protected:
			~residue_name_aligner_test_suite_fixture() noexcept = default;

		public:
			/// Check that residue_name_aligner::residue_name_align() gives the correct results for all permutations of the residues lists
			void check_residue_name_aligner_results(const residue_name_vec_vec &arg_residue_lists,
			                                        const bool_deq_vec         &arg_correct_presence_lists,
			                                        const size_vec_vec         &arg_correct_answer_lists
			                                        ) {
				do_check_residue_name_aligner( arg_residue_lists, arg_correct_presence_lists, arg_correct_answer_lists, false );
			}

			/// Check that residue_name_aligner::residue_name_align() throws an errro for all permutations of the residues lists
			void check_residue_name_aligner_throws(const residue_name_vec_vec &arg_residue_lists
			                                       ) {
				do_check_residue_name_aligner(arg_residue_lists, bool_deq_vec(), size_vec_vec(), true);
			}
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(residue_name_aligner_test_suite, cath::test::residue_name_aligner_test_suite_fixture)

// * Issues to test for multiples (more than two):
//    * do A-B, B-C, C-A, D [throw error]
//    * do A-B, C-D         [throw error]
//    * do A-B, B-C

/// \brief Check that the residue_name_aligner works correctly if two lists are completely identical
BOOST_AUTO_TEST_CASE(identical) {
	check_residue_name_aligner_results(
		{
			{ residue_name( 0 ), residue_name( 1 ), residue_name( 2 ), residue_name( 3 ) },
			{ residue_name( 0 ), residue_name( 1 ), residue_name( 2 ), residue_name( 3 ) }
		},
		{
			{ true, true, true, true },
			{ true, true, true, true }
		},
		{
			{ 0, 1, 2, 3 },
			{ 0, 1, 2, 3 }
		}
	);
}

/// \brief Check that the residue_name_aligner works correctly if the two lists each partially overlap with the other
BOOST_AUTO_TEST_CASE(partial_overlap) {
	check_residue_name_aligner_results(
		{
			{                    residue_name( 1 ), residue_name( 2 ), residue_name( 3 ) },
			{ residue_name( 0 ), residue_name( 1 ), residue_name( 2 )                    }
		},
		{
			{ false, true, true, true  },
			{ true,  true, true, false }
		},
		{
			{ 9999, 0, 1,    2 },
			{    0, 1, 2, 9999 }
		}
	);
}

/// \brief Check that the residue_name_aligner works correctly if one list is a (strict) subset of the other
BOOST_AUTO_TEST_CASE(strict_subset) {
	check_residue_name_aligner_results(
		{
			{                    residue_name( 1 ), residue_name( 2 ),                   },
			{ residue_name( 0 ), residue_name( 1 ), residue_name( 2 ), residue_name( 3 ) }
		},
		{
			{ false, true, true, false },
			{ true,  true, true, true  }
		},
		{
			{ 9999, 0, 1, 9999 },
			{    0, 1, 2,    3 }
		}
	);
}

/// \brief Check that the residue_name_aligner works correctly if the two lists are non consecutive (with gaps)
BOOST_AUTO_TEST_CASE(non_consecutive) {
	check_residue_name_aligner_results(
		{
			{ residue_name( 0 ), residue_name( 3 ),                    residue_name( 2 ), residue_name( 4 ), residue_name( 5 ) },
			{ residue_name( 0 ),                    residue_name( 1 ), residue_name( 2 ),                    residue_name( 5 ) }
		},
		{
			{ true,  true, false,  true,  true,  true },
			{ true, false,  true,  true, false,  true }
		},
		{
			{ 0,    1, 9999, 2, 3, 4 },
			{ 0, 9999,    1, 2, 0, 3 }
		}
	);
}

/// \brief Check that the residue_name_aligner correctly throws an error if the entries don't share any residue
BOOST_AUTO_TEST_CASE(out_of_order) {
	check_residue_name_aligner_throws(
		{
			{ residue_name( 2 ), residue_name( 1 ), residue_name( 0 ), residue_name( 3 ) },
			{ residue_name( 3 ), residue_name( 1 ), residue_name( 0 ), residue_name( 2 ) }
		}
	);
}

/// \brief Check that the residue_name_aligner correctly throws an error if the entries don't share any residues
///
/// I originally made this throw but after reflection, I couldn't see why this part should throw.
///
/// However, I do think that an attempt to superpose with fewer than three common coordinates should throw.
BOOST_AUTO_TEST_CASE(non_overlapping) {
	check_residue_name_aligner_results(
		{
			{ residue_name( 1 ), residue_name( 2 ), residue_name( 3 ),                                      },
			{                                                          residue_name( 4 ), residue_name( 0 ) }
		},
		{
			{ true,  true,  true,  false, false },
			{ false, false, false, true,  true  }
		},
		{
			{    0,    1,    2, 9999, 9999 },
			{ 9999, 9999, 9999,    0,    1 }
		}
	);
}

BOOST_AUTO_TEST_SUITE_END()
