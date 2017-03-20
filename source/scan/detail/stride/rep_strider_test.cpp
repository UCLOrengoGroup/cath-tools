/// \file
/// \brief The rep_strider test suite

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

#include "scan/detail/stride/rep_strider.hpp"

// #include "test/global_test_constants.hpp"

using namespace cath::scan::detail;
//using namespace std;

namespace cath {
	namespace test {

		/// \brief The rep_strider_test_suite_fixture to assist in testing rep_strider
		struct rep_strider_test_suite_fixture {
		protected:
			~rep_strider_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief Test suite to test the functionality of rep_strider
BOOST_FIXTURE_TEST_SUITE(rep_strider_test_suite, cath::test::rep_strider_test_suite_fixture)

/// \brief Test the ctor and getter with static_asserts
BOOST_AUTO_TEST_CASE(ctor_and_getter) {
	BOOST_CHECK( true );
	static_assert( rep_strider(   ).get_stride() == 0, "Check the ctor and getter work for default construction" );
	static_assert( rep_strider( 0 ).get_stride() == 0, "Check the ctor and getter work for construction from 0" );
	static_assert( rep_strider( 1 ).get_stride() == 1, "Check the ctor and getter work for construction from 1" );
	static_assert( rep_strider( 2 ).get_stride() == 2, "Check the ctor and getter work for construction from 2" );
	static_assert( rep_strider( 3 ).get_stride() == 3, "Check the ctor and getter work for construction from 3" );
	static_assert( rep_strider( 4 ).get_stride() == 4, "Check the ctor and getter work for construction from 4" );
}

/// \brief Test get_num_reps_of_num_residues_works() with static_asserts
BOOST_AUTO_TEST_CASE(get_num_reps_of_num_residues_works) {
	BOOST_CHECK( true );
	static_assert( get_num_reps_of_num_residues( rep_strider( 0 ), 0 ) == 0, "Check get_num_reps_of_num_residues_works() works for stride 0 and 0 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 0 ), 1 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 0 and 1 residues" );
	static_assert( get_num_reps_of_num_residues( rep_strider( 0 ), 2 ) == 2, "Check get_num_reps_of_num_residues_works() works for stride 0 and 2 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 0 ), 3 ) == 3, "Check get_num_reps_of_num_residues_works() works for stride 0 and 3 residue"  );

	static_assert( get_num_reps_of_num_residues( rep_strider( 1 ), 0 ) == 0, "Check get_num_reps_of_num_residues_works() works for stride 1 and 0 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 1 ), 1 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 1 and 1 residues" );
	static_assert( get_num_reps_of_num_residues( rep_strider( 1 ), 2 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 1 and 2 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 1 ), 3 ) == 2, "Check get_num_reps_of_num_residues_works() works for stride 1 and 3 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 1 ), 4 ) == 2, "Check get_num_reps_of_num_residues_works() works for stride 1 and 3 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 1 ), 5 ) == 3, "Check get_num_reps_of_num_residues_works() works for stride 1 and 3 residue"  );

	static_assert( get_num_reps_of_num_residues( rep_strider( 2 ), 0 ) == 0, "Check get_num_reps_of_num_residues_works() works for stride 2 and 0 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 2 ), 1 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 2 and 1 residues" );
	static_assert( get_num_reps_of_num_residues( rep_strider( 2 ), 2 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 2 and 2 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 2 ), 3 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 2 and 3 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 2 ), 4 ) == 2, "Check get_num_reps_of_num_residues_works() works for stride 2 and 3 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 2 ), 5 ) == 2, "Check get_num_reps_of_num_residues_works() works for stride 2 and 3 residue"  );

	static_assert( get_num_reps_of_num_residues( rep_strider( 3 ), 0 ) == 0, "Check get_num_reps_of_num_residues_works() works for stride 3 and 0 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 3 ), 1 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 3 and 1 residues" );
	static_assert( get_num_reps_of_num_residues( rep_strider( 3 ), 2 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 3 and 2 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 3 ), 3 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 3 and 3 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 3 ), 4 ) == 1, "Check get_num_reps_of_num_residues_works() works for stride 3 and 3 residue"  );
	static_assert( get_num_reps_of_num_residues( rep_strider( 3 ), 5 ) == 2, "Check get_num_reps_of_num_residues_works() works for stride 3 and 3 residue"  );
}

/// \brief Test get_index_of_rep_index() with static_asserts
BOOST_AUTO_TEST_CASE(get_index_of_rep_index_works) {
	BOOST_CHECK( true );

	static_assert( get_index_of_rep_index( rep_strider( 0 ), 0 ) ==  0, "Check get_index_of_rep_index() works for stride 0 and rep 0" );
	static_assert( get_index_of_rep_index( rep_strider( 0 ), 1 ) ==  1, "Check get_index_of_rep_index() works for stride 0 and rep 1" );
	static_assert( get_index_of_rep_index( rep_strider( 0 ), 2 ) ==  2, "Check get_index_of_rep_index() works for stride 0 and rep 2" );
	static_assert( get_index_of_rep_index( rep_strider( 0 ), 3 ) ==  3, "Check get_index_of_rep_index() works for stride 0 and rep 3" );

	static_assert( get_index_of_rep_index( rep_strider( 1 ), 0 ) ==  0, "Check get_index_of_rep_index() works for stride 1 and rep 0" );
	static_assert( get_index_of_rep_index( rep_strider( 1 ), 1 ) ==  2, "Check get_index_of_rep_index() works for stride 1 and rep 1" );
	static_assert( get_index_of_rep_index( rep_strider( 1 ), 2 ) ==  4, "Check get_index_of_rep_index() works for stride 1 and rep 2" );
	static_assert( get_index_of_rep_index( rep_strider( 1 ), 3 ) ==  6, "Check get_index_of_rep_index() works for stride 1 and rep 3" );
	static_assert( get_index_of_rep_index( rep_strider( 1 ), 4 ) ==  8, "Check get_index_of_rep_index() works for stride 1 and rep 3" );
	static_assert( get_index_of_rep_index( rep_strider( 1 ), 5 ) == 10, "Check get_index_of_rep_index() works for stride 1 and rep 3" );

	static_assert( get_index_of_rep_index( rep_strider( 2 ), 0 ) ==  0, "Check get_index_of_rep_index() works for stride 2 and rep 0" );
	static_assert( get_index_of_rep_index( rep_strider( 2 ), 1 ) ==  3, "Check get_index_of_rep_index() works for stride 2 and rep 1" );
	static_assert( get_index_of_rep_index( rep_strider( 2 ), 2 ) ==  6, "Check get_index_of_rep_index() works for stride 2 and rep 2" );
	static_assert( get_index_of_rep_index( rep_strider( 2 ), 3 ) ==  9, "Check get_index_of_rep_index() works for stride 2 and rep 3" );
	static_assert( get_index_of_rep_index( rep_strider( 2 ), 4 ) == 12, "Check get_index_of_rep_index() works for stride 2 and rep 3" );
	static_assert( get_index_of_rep_index( rep_strider( 2 ), 5 ) == 15, "Check get_index_of_rep_index() works for stride 2 and rep 3" );
 
 	static_assert( get_index_of_rep_index( rep_strider( 3 ), 0 ) ==  0, "Check get_index_of_rep_index() works for stride 3 and rep 0" );
 	static_assert( get_index_of_rep_index( rep_strider( 3 ), 1 ) ==  4, "Check get_index_of_rep_index() works for stride 3 and rep 1" );
 	static_assert( get_index_of_rep_index( rep_strider( 3 ), 2 ) ==  8, "Check get_index_of_rep_index() works for stride 3 and rep 2" );
 	static_assert( get_index_of_rep_index( rep_strider( 3 ), 3 ) == 12, "Check get_index_of_rep_index() works for stride 3 and rep 3" );
 	static_assert( get_index_of_rep_index( rep_strider( 3 ), 4 ) == 16, "Check get_index_of_rep_index() works for stride 3 and rep 3" );
	static_assert( get_index_of_rep_index( rep_strider( 3 ), 5 ) == 20, "Check get_index_of_rep_index() works for stride 3 and rep 3" );
}

//add_static_tests_here()
//
//inline constexpr rep_strider get_query_from_strider(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
//inline constexpr rep_strider get_query_to_strider(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
//inline constexpr rep_strider get_index_from_strider(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
//inline constexpr rep_strider get_index_to_strider(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
//inline constexpr index_type get_from_co_stride(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values
//inline constexpr index_type get_to_co_stride(const scan_stride &arg_stride ///< The scan_stride object containing the relevant stride values

BOOST_AUTO_TEST_SUITE_END()
