/// \file
/// \brief The cath_id_score_category test suite

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
#include <boost/utility/string_ref.hpp>

#include "resolve_hits/file/cath_id_score_category.hpp"

using namespace cath::rslv;

using boost::string_ref;
using std::string;

BOOST_AUTO_TEST_SUITE(cath_id_score_category_test_suite)

BOOST_AUTO_TEST_CASE(ending_round_1_is_norm) {
	const string id = "bob_round_1";
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, true  ), cath_id_score_category::NORMAL  );
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, false ), cath_id_score_category::NORMAL  );
}

BOOST_AUTO_TEST_CASE(ending_round_2_is_norm) {
	const string id = "bob_round_2";
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, true  ), cath_id_score_category::NORMAL  );
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, false ), cath_id_score_category::NORMAL  );
}

BOOST_AUTO_TEST_CASE(ending_round_74_is_norm) {
	const string id = "bob_round_74";
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, true  ), cath_id_score_category::NORMAL  );
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, false ), cath_id_score_category::NORMAL  );
}

BOOST_AUTO_TEST_CASE(id_1u5mA01_round_3_is_norm ) {
	const string id = "1u5mA01_round_3";
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, true  ), cath_id_score_category::NORMAL  );
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, false ), cath_id_score_category::NORMAL  );
}

BOOST_AUTO_TEST_CASE(dc_example_is_dc_type) {
	const string id = "dc_72a964d791dea7a3dd35a8bbf49385b8";
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, true  ), cath_id_score_category::DC_TYPE );
	BOOST_CHECK_EQUAL( cath_score_category_of_id( { id }, false ), cath_id_score_category::NORMAL  );
}

BOOST_AUTO_TEST_CASE(hmmer_evalues_are_suspicious_works) {
	static_assert( ! hmmer_evalues_are_suspicious( 0.0009, 0.0009 ), "" );
	static_assert( ! hmmer_evalues_are_suspicious( 0.0009, 0.0010 ), "" );
	static_assert(   hmmer_evalues_are_suspicious( 0.0009, 0.0011 ), "" );

	static_assert( ! hmmer_evalues_are_suspicious( 0.0010, 0.0009 ), "" );
	static_assert( ! hmmer_evalues_are_suspicious( 0.0010, 0.0010 ), "" );
	static_assert(   hmmer_evalues_are_suspicious( 0.0010, 0.0011 ), "" );

	static_assert( ! hmmer_evalues_are_suspicious( 0.0011, 0.0009 ), "" );
	static_assert( ! hmmer_evalues_are_suspicious( 0.0011, 0.0010 ), "" );
	static_assert( ! hmmer_evalues_are_suspicious( 0.0011, 0.0011 ), "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE(bitscore_divisor_works) {
	static_assert( bitscore_divisor( false, false ) == 1.0, "" );
	static_assert( bitscore_divisor( false, true  ) == 1.0, "" );
	static_assert( bitscore_divisor( true,  false ) == 1.0, "" );
	static_assert( bitscore_divisor( true,  true  ) == 4.0, "" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
