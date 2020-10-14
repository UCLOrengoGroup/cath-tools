/// \file
/// \brief The crh_single_output_spec test suite

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

#include "cath/resolve_hits/options/spec/crh_single_output_spec.hpp"

using namespace cath::rslv;

BOOST_AUTO_TEST_SUITE(crh_single_output_spec_test_suite)

BOOST_AUTO_TEST_CASE(get_deprecated_suggestion_works) {
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec()                                                          ), ""                                  );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec()                        .set_summarise           ( true ) ), "--summarise-to-file -"             );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec()                        .set_generate_html_output( true ) ), "--html-output-to-file -"           );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec()                        .set_json_output         ( true ) ), "--json-output-to-file -"           );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec().set_output_file( "abc")                                  ), "--quiet --hits-text-to-file abc"   );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec().set_output_file( "abc").set_summarise           ( true ) ), "--quiet --summarise-to-file abc"   );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec().set_output_file( "abc").set_generate_html_output( true ) ), "--quiet --html-output-to-file abc" );
	BOOST_CHECK_EQUAL( get_deprecated_suggestion_str( crh_single_output_spec().set_output_file( "abc").set_json_output         ( true ) ), "--quiet --json-output-to-file abc" );
}

BOOST_AUTO_TEST_SUITE_END()

