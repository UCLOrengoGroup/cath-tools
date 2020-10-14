/// \file
/// \brief The co_stride test suite

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

#include <boost/test/unit_test.hpp>

#include "cath/common/size_t_literal.hpp"
#include "cath/scan/detail/stride/co_stride.hpp"
#include "cath/test/global_test_constants.hpp"

using namespace ::cath::common;
using namespace ::cath::scan::detail;
using namespace ::cath::scan::detail::detail;
using namespace ::std;

namespace cath {
	namespace test {

		/// \brief The _test_suite_fixture to assist in testing 
		struct co_stride_test_suite_fixture {
		protected:
			~co_stride_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(co_stride_test_suite, cath::test::co_stride_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(co_stride_works) {
	static_assert( co_stride( 0_z, 0_z ) ==  0, "co_stride( 0, 0 ) should be equal to  0" );
	static_assert( co_stride( 0_z, 1_z ) ==  1, "co_stride( 0, 1 ) should be equal to  1" );
	static_assert( co_stride( 0_z, 2_z ) ==  2, "co_stride( 0, 2 ) should be equal to  2" );
	static_assert( co_stride( 0_z, 3_z ) ==  3, "co_stride( 0, 3 ) should be equal to  3" );
	static_assert( co_stride( 0_z, 4_z ) ==  4, "co_stride( 0, 4 ) should be equal to  4" );
	static_assert( co_stride( 0_z, 5_z ) ==  5, "co_stride( 0, 5 ) should be equal to  5" );
	static_assert( co_stride( 1_z, 0_z ) ==  1, "co_stride( 1, 0 ) should be equal to  1" );
	static_assert( co_stride( 1_z, 1_z ) ==  1, "co_stride( 1, 1 ) should be equal to  1" );
	static_assert( co_stride( 1_z, 2_z ) ==  5, "co_stride( 1, 2 ) should be equal to  5" );
	static_assert( co_stride( 1_z, 3_z ) ==  3, "co_stride( 1, 3 ) should be equal to  3" );
	static_assert( co_stride( 1_z, 4_z ) ==  9, "co_stride( 1, 4 ) should be equal to  9" );
	static_assert( co_stride( 1_z, 5_z ) ==  5, "co_stride( 1, 5 ) should be equal to  5" );
	static_assert( co_stride( 2_z, 0_z ) ==  2, "co_stride( 2, 0 ) should be equal to  2" );
	static_assert( co_stride( 2_z, 1_z ) ==  5, "co_stride( 2, 1 ) should be equal to  5" );
	static_assert( co_stride( 2_z, 2_z ) ==  2, "co_stride( 2, 2 ) should be equal to  2" );
	static_assert( co_stride( 2_z, 3_z ) == 11, "co_stride( 2, 3 ) should be equal to 11" );
	static_assert( co_stride( 2_z, 4_z ) == 14, "co_stride( 2, 4 ) should be equal to 14" );
	static_assert( co_stride( 2_z, 5_z ) ==  5, "co_stride( 2, 5 ) should be equal to  5" );
	static_assert( co_stride( 3_z, 0_z ) ==  3, "co_stride( 3, 0 ) should be equal to  3" );
	static_assert( co_stride( 3_z, 1_z ) ==  3, "co_stride( 3, 1 ) should be equal to  3" );
	static_assert( co_stride( 3_z, 2_z ) == 11, "co_stride( 3, 2 ) should be equal to 11" );
	static_assert( co_stride( 3_z, 3_z ) ==  3, "co_stride( 3, 3 ) should be equal to  3" );
	static_assert( co_stride( 3_z, 4_z ) == 19, "co_stride( 3, 4 ) should be equal to 19" );
	static_assert( co_stride( 3_z, 5_z ) == 11, "co_stride( 3, 5 ) should be equal to 11" );
	static_assert( co_stride( 4_z, 0_z ) ==  4, "co_stride( 4, 0 ) should be equal to  4" );
	static_assert( co_stride( 4_z, 1_z ) ==  9, "co_stride( 4, 1 ) should be equal to  9" );
	static_assert( co_stride( 4_z, 2_z ) == 14, "co_stride( 4, 2 ) should be equal to 14" );
	static_assert( co_stride( 4_z, 3_z ) == 19, "co_stride( 4, 3 ) should be equal to 19" );
	static_assert( co_stride( 4_z, 4_z ) ==  4, "co_stride( 4, 4 ) should be equal to  4" );
	static_assert( co_stride( 4_z, 5_z ) == 29, "co_stride( 4, 5 ) should be equal to 29" );
	static_assert( co_stride( 5_z, 0_z ) ==  5, "co_stride( 5, 0 ) should be equal to  5" );
	static_assert( co_stride( 5_z, 1_z ) ==  5, "co_stride( 5, 1 ) should be equal to  5" );
	static_assert( co_stride( 5_z, 2_z ) ==  5, "co_stride( 5, 2 ) should be equal to  5" );
	static_assert( co_stride( 5_z, 3_z ) == 11, "co_stride( 5, 3 ) should be equal to 11" );
	static_assert( co_stride( 5_z, 4_z ) == 29, "co_stride( 5, 4 ) should be equal to 29" );
	static_assert( co_stride( 5_z, 5_z ) ==  5, "co_stride( 5, 5 ) should be equal to  5" );
	BOOST_CHECK( true );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(num_stride_neighbour_range_works) {
	static_assert( num_in_stride_neighbour_range( 0_z ) ==  1, "num_in_stride_neighbour_range( 0 ) should be equal to  1" );
	static_assert( num_in_stride_neighbour_range( 1_z ) ==  2, "num_in_stride_neighbour_range( 1 ) should be equal to  2" );
	static_assert( num_in_stride_neighbour_range( 2_z ) ==  3, "num_in_stride_neighbour_range( 2 ) should be equal to  3" );
	static_assert( num_in_stride_neighbour_range( 3_z ) ==  4, "num_in_stride_neighbour_range( 3 ) should be equal to  4" );
	static_assert( num_in_stride_neighbour_range( 4_z ) ==  5, "num_in_stride_neighbour_range( 4 ) should be equal to  5" );
	static_assert( num_in_stride_neighbour_range( 5_z ) ==  6, "num_in_stride_neighbour_range( 5 ) should be equal to  6" );
	static_assert( num_in_stride_neighbour_range( 6_z ) ==  7, "num_in_stride_neighbour_range( 6 ) should be equal to  7" );
	static_assert( num_in_stride_neighbour_range( 7_z ) ==  8, "num_in_stride_neighbour_range( 7 ) should be equal to  8" );
	static_assert( num_in_stride_neighbour_range( 8_z ) ==  9, "num_in_stride_neighbour_range( 8 ) should be equal to  9" );
	static_assert( num_in_stride_neighbour_range( 9_z ) == 10, "num_in_stride_neighbour_range( 9 ) should be equal to 10" );
	BOOST_CHECK( true );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(stride_neighbours_index_of_centre_works) {
	static_assert( stride_neighbour_index_of_centre( 0_z ) == 0, "stride_neighbour_index_of_centre( 0 ) should be equal to 0" );
	static_assert( stride_neighbour_index_of_centre( 1_z ) == 0, "stride_neighbour_index_of_centre( 1 ) should be equal to 0" );
	static_assert( stride_neighbour_index_of_centre( 2_z ) == 1, "stride_neighbour_index_of_centre( 2 ) should be equal to 1" );
	static_assert( stride_neighbour_index_of_centre( 3_z ) == 1, "stride_neighbour_index_of_centre( 3 ) should be equal to 1" );
	static_assert( stride_neighbour_index_of_centre( 4_z ) == 2, "stride_neighbour_index_of_centre( 4 ) should be equal to 2" );
	static_assert( stride_neighbour_index_of_centre( 5_z ) == 2, "stride_neighbour_index_of_centre( 5 ) should be equal to 2" );
	static_assert( stride_neighbour_index_of_centre( 6_z ) == 3, "stride_neighbour_index_of_centre( 6 ) should be equal to 3" );
	static_assert( stride_neighbour_index_of_centre( 7_z ) == 3, "stride_neighbour_index_of_centre( 7 ) should be equal to 3" );
	static_assert( stride_neighbour_index_of_centre( 8_z ) == 4, "stride_neighbour_index_of_centre( 8 ) should be equal to 4" );
	static_assert( stride_neighbour_index_of_centre( 9_z ) == 4, "stride_neighbour_index_of_centre( 9 ) should be equal to 4" );
	BOOST_CHECK( true );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(entry_index_of_stride_neighbour_index_works) {
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  0_z,  0_z,  3_z ) == make_pair(  true,  0_z ), "entry_index_of_stride_neighbour_index( 0_z,  0_z,  0_z,  3_z ) should be equal to    0" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  0_z,  1_z,  3_z ) == make_pair(  true,  1_z ), "entry_index_of_stride_neighbour_index( 0_z,  0_z,  1_z,  3_z ) should be equal to    1" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  0_z,  2_z,  3_z ) == make_pair(  true,  2_z ), "entry_index_of_stride_neighbour_index( 0_z,  0_z,  2_z,  3_z ) should be equal to    2" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  0_z,  3_z,  3_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 0_z,  0_z,  3_z,  3_z ) should be equal to none" );

	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  1_z,  0_z,  5_z ) == make_pair(  true,  0_z ), "entry_index_of_stride_neighbour_index( 0_z,  1_z,  0_z,  5_z ) should be equal to    0" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  1_z,  0_z,  5_z ) == make_pair(  true,  1_z ), "entry_index_of_stride_neighbour_index( 1_z,  1_z,  0_z,  5_z ) should be equal to    1" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  1_z,  2_z,  5_z ) == make_pair(  true,  2_z ), "entry_index_of_stride_neighbour_index( 0_z,  1_z,  2_z,  5_z ) should be equal to    2" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  1_z,  2_z,  5_z ) == make_pair(  true,  3_z ), "entry_index_of_stride_neighbour_index( 1_z,  1_z,  2_z,  5_z ) should be equal to    3" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  1_z,  4_z,  5_z ) == make_pair(  true,  4_z ), "entry_index_of_stride_neighbour_index( 0_z,  1_z,  4_z,  5_z ) should be equal to    4" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  1_z,  4_z,  5_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 1_z,  1_z,  4_z,  5_z ) should be equal to none" );

	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  2_z,  0_z,  7_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 0_z,  2_z,  0_z,  7_z ) should be equal to none" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  2_z,  0_z,  7_z ) == make_pair(  true,  0_z ), "entry_index_of_stride_neighbour_index( 1_z,  2_z,  0_z,  7_z ) should be equal to    0" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  2_z,  0_z,  7_z ) == make_pair(  true,  1_z ), "entry_index_of_stride_neighbour_index( 2_z,  2_z,  0_z,  7_z ) should be equal to    1" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  2_z,  3_z,  7_z ) == make_pair(  true,  2_z ), "entry_index_of_stride_neighbour_index( 0_z,  2_z,  3_z,  7_z ) should be equal to    2" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  2_z,  3_z,  7_z ) == make_pair(  true,  3_z ), "entry_index_of_stride_neighbour_index( 1_z,  2_z,  3_z,  7_z ) should be equal to    3" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  2_z,  3_z,  7_z ) == make_pair(  true,  4_z ), "entry_index_of_stride_neighbour_index( 2_z,  2_z,  3_z,  7_z ) should be equal to    4" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  2_z,  6_z,  7_z ) == make_pair(  true,  5_z ), "entry_index_of_stride_neighbour_index( 0_z,  2_z,  6_z,  7_z ) should be equal to    5" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  2_z,  6_z,  7_z ) == make_pair(  true,  6_z ), "entry_index_of_stride_neighbour_index( 1_z,  2_z,  6_z,  7_z ) should be equal to    6" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  2_z,  6_z,  7_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 2_z,  2_z,  6_z,  7_z ) should be equal to none" );

	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  3_z,  0_z,  9_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 0_z,  3_z,  0_z,  9_z ) should be equal to none" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  3_z,  0_z,  9_z ) == make_pair(  true,  0_z ), "entry_index_of_stride_neighbour_index( 1_z,  3_z,  0_z,  9_z ) should be equal to    0" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  3_z,  0_z,  9_z ) == make_pair(  true,  1_z ), "entry_index_of_stride_neighbour_index( 2_z,  3_z,  0_z,  9_z ) should be equal to    1" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 3_z,  3_z,  0_z,  9_z ) == make_pair(  true,  2_z ), "entry_index_of_stride_neighbour_index( 3_z,  3_z,  0_z,  9_z ) should be equal to    2" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  3_z,  4_z,  9_z ) == make_pair(  true,  3_z ), "entry_index_of_stride_neighbour_index( 0_z,  3_z,  4_z,  9_z ) should be equal to    3" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  3_z,  4_z,  9_z ) == make_pair(  true,  4_z ), "entry_index_of_stride_neighbour_index( 1_z,  3_z,  4_z,  9_z ) should be equal to    4" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  3_z,  4_z,  9_z ) == make_pair(  true,  5_z ), "entry_index_of_stride_neighbour_index( 2_z,  3_z,  4_z,  9_z ) should be equal to    5" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 3_z,  3_z,  4_z,  9_z ) == make_pair(  true,  6_z ), "entry_index_of_stride_neighbour_index( 3_z,  3_z,  4_z,  9_z ) should be equal to    6" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  3_z,  8_z,  9_z ) == make_pair(  true,  7_z ), "entry_index_of_stride_neighbour_index( 0_z,  3_z,  8_z,  9_z ) should be equal to    7" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  3_z,  8_z,  9_z ) == make_pair(  true,  8_z ), "entry_index_of_stride_neighbour_index( 1_z,  3_z,  8_z,  9_z ) should be equal to    8" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  3_z,  8_z,  9_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 2_z,  3_z,  8_z,  9_z ) should be equal to none" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 3_z,  3_z,  8_z,  9_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 3_z,  3_z,  8_z,  9_z ) should be equal to none" );

	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  4_z,  0_z, 11_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 0_z,  4_z,  0_z,  9_z ) should be equal to none" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  4_z,  0_z, 11_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 1_z,  4_z,  0_z,  9_z ) should be equal to none" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  4_z,  0_z, 11_z ) == make_pair(  true,  0_z ), "entry_index_of_stride_neighbour_index( 2_z,  4_z,  0_z,  9_z ) should be equal to    0" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 3_z,  4_z,  0_z, 11_z ) == make_pair(  true,  1_z ), "entry_index_of_stride_neighbour_index( 3_z,  4_z,  0_z,  9_z ) should be equal to    1" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 4_z,  4_z,  0_z, 11_z ) == make_pair(  true,  2_z ), "entry_index_of_stride_neighbour_index( 4_z,  4_z,  0_z,  9_z ) should be equal to    2" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  4_z,  5_z, 11_z ) == make_pair(  true,  3_z ), "entry_index_of_stride_neighbour_index( 0_z,  4_z,  5_z,  9_z ) should be equal to    3" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  4_z,  5_z, 11_z ) == make_pair(  true,  4_z ), "entry_index_of_stride_neighbour_index( 1_z,  4_z,  5_z,  9_z ) should be equal to    4" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  4_z,  5_z, 11_z ) == make_pair(  true,  5_z ), "entry_index_of_stride_neighbour_index( 2_z,  4_z,  5_z,  9_z ) should be equal to    5" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 3_z,  4_z,  5_z, 11_z ) == make_pair(  true,  6_z ), "entry_index_of_stride_neighbour_index( 3_z,  4_z,  5_z,  9_z ) should be equal to    6" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 4_z,  4_z,  5_z, 11_z ) == make_pair(  true,  7_z ), "entry_index_of_stride_neighbour_index( 4_z,  4_z,  5_z,  9_z ) should be equal to    7" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 0_z,  4_z, 10_z, 11_z ) == make_pair(  true,  8_z ), "entry_index_of_stride_neighbour_index( 0_z,  4_z, 10_z,  9_z ) should be equal to    8" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 1_z,  4_z, 10_z, 11_z ) == make_pair(  true,  9_z ), "entry_index_of_stride_neighbour_index( 1_z,  4_z, 10_z,  9_z ) should be equal to    9" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 2_z,  4_z, 10_z, 11_z ) == make_pair(  true, 10_z ), "entry_index_of_stride_neighbour_index( 2_z,  4_z, 10_z,  9_z ) should be equal to   10" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 3_z,  4_z, 10_z, 11_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 3_z,  4_z, 10_z,  9_z ) should be equal to none" );
	static_assert( entry_index_of_stride_neighbour_index_impl( 4_z,  4_z, 10_z, 11_z ) == make_pair( false,  0_z ), "entry_index_of_stride_neighbour_index( 4_z,  4_z, 10_z,  9_z ) should be equal to none" );
	BOOST_CHECK( true );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(entry_index_of_stride_rep_works) {
	static_assert( entry_index_of_stride_rep(  0_z, 0_z ) ==  0_z, "entry_index_of_stride_rep(  0_z, 0_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  1_z, 0_z ) ==  1_z, "entry_index_of_stride_rep(  1_z, 0_z ) should be equal to  1" );
	static_assert( entry_index_of_stride_rep(  2_z, 0_z ) ==  2_z, "entry_index_of_stride_rep(  2_z, 0_z ) should be equal to  2" );

	static_assert( entry_index_of_stride_rep(  0_z, 1_z ) ==  0_z, "entry_index_of_stride_rep(  0_z, 1_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  1_z, 1_z ) ==  0_z, "entry_index_of_stride_rep(  1_z, 1_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  2_z, 1_z ) ==  2_z, "entry_index_of_stride_rep(  2_z, 1_z ) should be equal to  2" );
	static_assert( entry_index_of_stride_rep(  3_z, 1_z ) ==  2_z, "entry_index_of_stride_rep(  3_z, 1_z ) should be equal to  2" );
	static_assert( entry_index_of_stride_rep(  4_z, 1_z ) ==  4_z, "entry_index_of_stride_rep(  4_z, 1_z ) should be equal to  4" );

	static_assert( entry_index_of_stride_rep(  0_z, 2_z ) ==  0_z, "entry_index_of_stride_rep(  0_z, 2_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  1_z, 2_z ) ==  0_z, "entry_index_of_stride_rep(  1_z, 2_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  2_z, 2_z ) ==  3_z, "entry_index_of_stride_rep(  2_z, 2_z ) should be equal to  3" );
	static_assert( entry_index_of_stride_rep(  3_z, 2_z ) ==  3_z, "entry_index_of_stride_rep(  3_z, 2_z ) should be equal to  3" );
	static_assert( entry_index_of_stride_rep(  4_z, 2_z ) ==  3_z, "entry_index_of_stride_rep(  4_z, 2_z ) should be equal to  3" );
	static_assert( entry_index_of_stride_rep(  5_z, 2_z ) ==  6_z, "entry_index_of_stride_rep(  5_z, 2_z ) should be equal to  6" );
	static_assert( entry_index_of_stride_rep(  6_z, 2_z ) ==  6_z, "entry_index_of_stride_rep(  6_z, 2_z ) should be equal to  6" );

	static_assert( entry_index_of_stride_rep(  0_z, 3_z ) ==  0_z, "entry_index_of_stride_rep(  0_z, 3_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  1_z, 3_z ) ==  0_z, "entry_index_of_stride_rep(  1_z, 3_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  2_z, 3_z ) ==  0_z, "entry_index_of_stride_rep(  2_z, 3_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  3_z, 3_z ) ==  4_z, "entry_index_of_stride_rep(  3_z, 3_z ) should be equal to  4" );
	static_assert( entry_index_of_stride_rep(  4_z, 3_z ) ==  4_z, "entry_index_of_stride_rep(  4_z, 3_z ) should be equal to  4" );
	static_assert( entry_index_of_stride_rep(  5_z, 3_z ) ==  4_z, "entry_index_of_stride_rep(  5_z, 3_z ) should be equal to  4" );
	static_assert( entry_index_of_stride_rep(  6_z, 3_z ) ==  4_z, "entry_index_of_stride_rep(  6_z, 3_z ) should be equal to  4" );
	static_assert( entry_index_of_stride_rep(  7_z, 3_z ) ==  8_z, "entry_index_of_stride_rep(  7_z, 3_z ) should be equal to  8" );
	static_assert( entry_index_of_stride_rep(  8_z, 3_z ) ==  8_z, "entry_index_of_stride_rep(  8_z, 3_z ) should be equal to  8" );

	static_assert( entry_index_of_stride_rep(  0_z, 4_z ) ==  0_z, "entry_index_of_stride_rep(  0_z, 4_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  1_z, 4_z ) ==  0_z, "entry_index_of_stride_rep(  1_z, 4_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  2_z, 4_z ) ==  0_z, "entry_index_of_stride_rep(  2_z, 4_z ) should be equal to  0" );
	static_assert( entry_index_of_stride_rep(  3_z, 4_z ) ==  5_z, "entry_index_of_stride_rep(  3_z, 4_z ) should be equal to  5" );
	static_assert( entry_index_of_stride_rep(  4_z, 4_z ) ==  5_z, "entry_index_of_stride_rep(  4_z, 4_z ) should be equal to  5" );
	static_assert( entry_index_of_stride_rep(  5_z, 4_z ) ==  5_z, "entry_index_of_stride_rep(  5_z, 4_z ) should be equal to  5" );
	static_assert( entry_index_of_stride_rep(  6_z, 4_z ) ==  5_z, "entry_index_of_stride_rep(  6_z, 4_z ) should be equal to  5" );
	static_assert( entry_index_of_stride_rep(  7_z, 4_z ) ==  5_z, "entry_index_of_stride_rep(  7_z, 4_z ) should be equal to  5" );
	static_assert( entry_index_of_stride_rep(  8_z, 4_z ) == 10_z, "entry_index_of_stride_rep(  8_z, 4_z ) should be equal to 10" );
	static_assert( entry_index_of_stride_rep(  9_z, 4_z ) == 10_z, "entry_index_of_stride_rep(  9_z, 4_z ) should be equal to 10" );
	static_assert( entry_index_of_stride_rep( 10_z, 4_z ) == 10_z, "entry_index_of_stride_rep( 10_z, 4_z ) should be equal to 10" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
