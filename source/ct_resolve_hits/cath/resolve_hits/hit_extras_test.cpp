/// \file
/// \brief The hit_extras test suite

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

#include <boost/optional/optional_io.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/common/type_aliases.hpp"
#include "cath/resolve_hits/hit_extras.hpp"

using namespace cath;
using namespace cath::rslv;

using boost::make_optional;
using boost::none;

BOOST_AUTO_TEST_SUITE(hit_extras_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	hit_extras_store the_store;
	
	the_store.push_back< hit_extra_cat::COND_EVAL >( 4.0 );
	the_store.push_back< hit_extra_cat::ALND_RGNS >( "#" );
	the_store.push_back< hit_extra_cat::INDP_EVAL >( 8.0 );

	BOOST_CHECK_EQUAL( to_string( the_store ), R"(hit_extras_store[cond-evalue:4.000000,aligned-regions:#,indp-evalue:8.000000])" );
}

BOOST_AUTO_TEST_CASE(get_first_works) {
	hit_extras_store the_store;
	the_store.push_back< hit_extra_cat::COND_EVAL >( 4.0 );
	BOOST_CHECK_EQUAL( get_first<hit_extra_cat::COND_EVAL>( the_store ), make_optional( 4.0 ) );
	BOOST_CHECK_EQUAL( get_first<hit_extra_cat::INDP_EVAL>( the_store ), doub_opt{ none }     );
}

BOOST_AUTO_TEST_SUITE_END()
