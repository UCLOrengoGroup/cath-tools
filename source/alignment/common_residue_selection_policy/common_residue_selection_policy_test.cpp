/// \file
/// \brief The common_residue_selection_policy test suite

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

#include <boost/mpl/list.hpp>

#include "alignment/alignment.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_min_score_policy.h"
#include "alignment/pair_alignment.h"
#include "common/boost_addenda/test/boost_check_equal_ranges.h"
#include "common/boost_check_no_throw_diag.h"
#include "common/size_t_literal.h"
#include "common/test_tools.h"
#include "exception/invalid_argument_exception.h"
#include "test/global_test_constants.h"

#include <limits>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::common::test;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The common_residue_selection_policy_test_suite_fixture to assist in testing common_residue_selection_policy
		struct common_residue_selection_policy_test_suite_fixture : public global_test_constants {
		protected:
			~common_residue_selection_policy_test_suite_fixture() noexcept = default;

		public:
			void check_policy_on_score_aln(const common_residue_selection_policy &,
			                               const size_vec &) const;

			const opt_aln_posn_vec aln_list_a       = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
			                                            11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

			const opt_aln_posn_vec aln_list_b       = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
			                                            11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

			const alignment        unscored_aln     = { alignment_offset_1_factory( { aln_list_a, aln_list_b } ) };

			const opt_score_vec    alignment_scores = { 47.7, 1, 94.2, 94.9, 52.8, 35.2, 67.3, 29, 12.8, 30.9,
			                                            27.6, 79.9, 35, 8.3, 21.6,88.1, 4.4, 45.7, 9.1, 95.9 };

			alignment              scored_aln       = { set_pair_alignment_duplicate_scores_copy(unscored_aln, alignment_scores) };
		};

	}
}

/// \brief TODOCUMENT
void cath::test::common_residue_selection_policy_test_suite_fixture::check_policy_on_score_aln(const common_residue_selection_policy &arg_policy,          ///< TODOCUMENT
                                                                                               const size_vec                        &arg_expected_indices ///< TODOCUMENT
                                                                                               ) const {
	const alignment &aln_const_ref(scored_aln);
	const size_vec   got_indices = select_common_residues_of_pair_alignment( arg_policy, aln_const_ref );
	BOOST_CHECK_EQUAL_RANGES( arg_expected_indices, got_indices );
}

BOOST_FIXTURE_TEST_SUITE(common_residue_selection_policy_test_suite, cath::test::common_residue_selection_policy_test_suite_fixture)

/// \brief A type-list containing all common_residue_selection_policy types that don't use scores.
using all_common_coord_non_score_based_selection_policy_types = boost::mpl::list<common_residue_select_all_policy>;

/// \brief A type-list containing all common_residue_selection_policy types that do use scores.
using all_common_coord_score_based_selection_policy_types = boost::mpl::list<common_residue_select_best_score_percent_policy, common_residue_select_min_score_policy>;

/// \brief A type-list containing all common_residue_selection_policy types.
using all_common_residue_selection_policy_types = boost::mpl::list<common_residue_select_all_policy, common_residue_select_best_score_percent_policy, common_residue_select_min_score_policy>;

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(constructor, common_residue_selection_policy_type, all_common_residue_selection_policy_types) {
	const common_residue_selection_policy_type my_policy{};
	BOOST_CHECK_NO_THROW_DIAG(const common_residue_selection_policy_type my_policy{});
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(does_not_throw_for_unscored_alignment, common_residue_selection_policy_type, all_common_coord_non_score_based_selection_policy_types) {
	const common_residue_selection_policy_type my_policy{};
	const alignment &aln_const_ref(unscored_aln);
	BOOST_CHECK_NO_THROW_DIAG(const vector<alignment::size_type> common_residue_indices = my_policy.select_common_residues(aln_const_ref, 0, 1));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(throws_for_unscored_alignment, common_residue_selection_policy_type, all_common_coord_score_based_selection_policy_types) {
	const common_residue_selection_policy_type my_policy{};
	const alignment &aln_const_ref(unscored_aln);
	BOOST_CHECK_THROW(const vector<alignment::size_type> common_residue_indices = my_policy.select_common_residues(aln_const_ref, 0, 1), invalid_argument_exception);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(throws_for_entry_out_of_range, common_residue_selection_policy_type, all_common_residue_selection_policy_types) {
	const common_residue_selection_policy_type my_policy{};
	const alignment &aln_const_ref(scored_aln);
	BOOST_CHECK_THROW(const vector<alignment::size_type> common_residue_indices = my_policy.select_common_residues(aln_const_ref, 0, 2), invalid_argument_exception);
	BOOST_CHECK_THROW(const vector<alignment::size_type> common_residue_indices = my_policy.select_common_residues(aln_const_ref, 2, 0), invalid_argument_exception);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(does_not_throw_for_scored_alignment, common_residue_selection_policy_type, all_common_residue_selection_policy_types) {
	const common_residue_selection_policy_type my_policy{};
	const alignment &aln_const_ref(scored_aln);
	BOOST_CHECK_NO_THROW_DIAG(const vector<alignment::size_type> common_residue_indices = my_policy.select_common_residues(aln_const_ref, 0, 1));
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(handles_an_empty_alignment, common_residue_selection_policy_type, all_common_residue_selection_policy_types) {
	alignment empty_pair_alignment(alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT);
	set_pair_alignment_duplicate_scores( empty_pair_alignment, opt_score_vec() );
	const common_residue_selection_policy_type my_policy{};
	const alignment &aln_const_ref(empty_pair_alignment);
	const vector<alignment::size_type> common_residue_indices = select_common_residues_of_pair_alignment(my_policy, aln_const_ref);
	BOOST_CHECK_EQUAL(common_residue_indices.size(), 0_z);
}

BOOST_AUTO_TEST_CASE(best_score_percent_ctor_throws_with_invalid_args) {
	BOOST_CHECK_THROW( common_residue_select_best_score_percent_policy pol(  -0.001                ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_best_score_percent_policy pol( 100.001                ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_best_score_percent_policy pol( DOUBLE_INFINITY()      ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_best_score_percent_policy pol( DOUBLE_QUIET_NAN ()    ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_best_score_percent_policy pol( DOUBLE_SIGNALING_NAN() ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_best_score_percent_policy pol( DOUBLE_MAX()           ), invalid_argument_exception );
}

BOOST_AUTO_TEST_CASE(min_score_ctor_throws_with_invalid_args) {
	BOOST_CHECK_THROW( common_residue_select_min_score_policy pol(  -0.5                  ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_min_score_policy pol( 100.001                ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_min_score_policy pol( DOUBLE_INFINITY()      ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_min_score_policy pol( DOUBLE_QUIET_NAN()     ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_min_score_policy pol( DOUBLE_SIGNALING_NAN() ), invalid_argument_exception );
	BOOST_CHECK_THROW( common_residue_select_min_score_policy pol( DOUBLE_MAX()           ), invalid_argument_exception );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(all) {
	check_policy_on_score_aln(
		common_residue_select_all_policy(),
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(best_score_percent_0) {
	check_policy_on_score_aln(
		common_residue_select_best_score_percent_policy(0),
		vector<alignment::size_type>()
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(best_score_percent_90) {
	check_policy_on_score_aln(
		common_residue_select_best_score_percent_policy(90),
		{ 0, 2, 3, 4, 5, 6, 9, 11, 12, 15, 17, 19 }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(best_score_percent_default) {
	check_policy_on_score_aln(
		common_residue_select_best_score_percent_policy(),
		{ 2, 3, 4, 6, 11, 15, 19 }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(best_score_percent_100) {
	check_policy_on_score_aln(
		common_residue_select_best_score_percent_policy(100.0),
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }
	);
}



/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(min_score_0) {
	check_policy_on_score_aln(
		common_residue_select_min_score_policy(0.0),
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(min_score_default) {
	check_policy_on_score_aln(
		common_residue_select_min_score_policy(),
		{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(min_score_35) {
	check_policy_on_score_aln(
		common_residue_select_min_score_policy(35.0),
		{ 0, 2, 3, 4, 5, 6, 11, 15, 17, 19 }
	);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(min_score_100) {
	check_policy_on_score_aln(
		common_residue_select_min_score_policy(100.0),
		vector<alignment::size_type>()
	);
}

/// \brief Check operator== and operator!= work over a range of different common_residue_selection_policies
BOOST_AUTO_TEST_CASE(equality_and_inequality_operators_work) {
	check_equality_operators_on_diff_vals_range( get_all_common_residue_selection_policies() );
}

/// \brief TODOCUMENT
///
/// \todo Improve the standard tools for testing less-than and then apply them here
BOOST_AUTO_TEST_CASE(less_than_operator_is_sensible) {
	const auto common_residue_selection_policies = get_all_common_residue_selection_policies();
	for (const auto &x : common_residue_selection_policies) {
		for (const auto &y : common_residue_selection_policies) {
			if ( x == y ) {
				BOOST_CHECK( !( x < y) );
				BOOST_CHECK( !( y < x) );
			}
			else {
				BOOST_CHECK(   ( x < y) ||   ( y < x) );
				BOOST_CHECK( ! ( x < y) || ! ( y < x) );
			}
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()
