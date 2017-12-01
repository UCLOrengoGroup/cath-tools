/// \file
/// \brief The align_refining test suite

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

#include "align_refining.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

using namespace cath::align;

using boost::lexical_cast;
using std::string;

BOOST_AUTO_TEST_SUITE(align_refining_test_suite)

BOOST_AUTO_TEST_CASE(to_string_works) {
	BOOST_TEST( to_string( align_refining::NO    ) == "NO"    );
	BOOST_TEST( to_string( align_refining::LIGHT ) == "LIGHT" );
	BOOST_TEST( to_string( align_refining::HEAVY ) == "HEAVY" );
}

BOOST_AUTO_TEST_CASE(insertion_operator_works) {
	BOOST_TEST( lexical_cast<string>( align_refining::NO    ) == "NO"    );
	BOOST_TEST( lexical_cast<string>( align_refining::LIGHT ) == "LIGHT" );
	BOOST_TEST( lexical_cast<string>( align_refining::HEAVY ) == "HEAVY" );
}

BOOST_AUTO_TEST_CASE(extraction_operator_works) {
	BOOST_TEST( lexical_cast<align_refining>( "NO"    ) == align_refining::NO    );
	BOOST_TEST( lexical_cast<align_refining>( "no"    ) == align_refining::NO    );

	BOOST_TEST( lexical_cast<align_refining>( "LIGHT" ) == align_refining::LIGHT );
	BOOST_TEST( lexical_cast<align_refining>( "light" ) == align_refining::LIGHT );

	BOOST_TEST( lexical_cast<align_refining>( "HEAVY" ) == align_refining::HEAVY );
	BOOST_TEST( lexical_cast<align_refining>( "heavy" ) == align_refining::HEAVY );
}

BOOST_AUTO_TEST_SUITE_END()
