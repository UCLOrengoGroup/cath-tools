/// \file
/// \brief The trim_spec test suite

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

#include "cath/resolve_hits/trim/trim_spec.hpp"

using namespace cath::rslv;
using namespace cath::seq;

BOOST_AUTO_TEST_SUITE(trim_spec_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert( trim_spec{ 5, 3 }.get_full_length   () == 5, "" );
	static_assert( trim_spec{ 5, 3 }.get_total_trimming() == 3, "" );

	// static_assert( trim_spec{ 2, 2 }.get_full_length() == 2, "" ); //Invalid

	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 0 ) == 0, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 1 ) == 0, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 2 ) == 0, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 3 ) == 1, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 4 ) == 2, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 5 ) == 3, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 6 ) == 3, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 7 ) == 3, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 8 ) == 3, "" );
	static_assert( total_trimming_of_length( trim_spec{ 5, 3 }, 9 ) == 3, "" );

	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 0 ) == 0, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 1 ) == 1, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 2 ) == 2, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 3 ) == 2, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 4 ) == 2, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 5 ) == 2, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 6 ) == 3, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 7 ) == 4, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 8 ) == 5, "" );
	static_assert( length_after_trimming   ( trim_spec{ 5, 3 }, 9 ) == 6, "" );

	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 0 ) == 0, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 1 ) == 0, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 2 ) == 0, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 3 ) == 0, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 4 ) == 1, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 5 ) == 1, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 6 ) == 1, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 7 ) == 1, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 8 ) == 1, "" );
	static_assert( start_trimming_of_length( trim_spec{ 5, 3 }, 9 ) == 1, "" );

	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 0 ) == 0, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 1 ) == 0, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 2 ) == 0, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 3 ) == 1, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 4 ) == 1, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 5 ) == 2, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 6 ) == 2, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 7 ) == 2, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 8 ) == 2, "" );
	static_assert( stop_trimming_of_length ( trim_spec{ 5, 3 }, 9 ) == 2, "" );

	// static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 10 ) == 0, "" ); // Invalid

	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 11 ) == residx_residx_pair{ 11, 11 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 12 ) == residx_residx_pair{ 11, 12 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 13 ) == residx_residx_pair{ 11, 12 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 14 ) == residx_residx_pair{ 12, 13 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 15 ) == residx_residx_pair{ 12, 13 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 16 ) == residx_residx_pair{ 12, 14 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 17 ) == residx_residx_pair{ 12, 15 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 18 ) == residx_residx_pair{ 12, 16 }, "" );
	static_assert( trim_copy_start_stop    ( trim_spec{ 5, 3 }, 11, 19 ) == residx_residx_pair{ 12, 17 }, "" );

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
