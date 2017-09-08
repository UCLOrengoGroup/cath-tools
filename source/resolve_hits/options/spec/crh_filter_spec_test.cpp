/// \file
/// \brief The crh_filter_spec test suite

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
#include <boost/test/auto_unit_test.hpp>

#include "resolve_hits/options/spec/crh_filter_spec.hpp"

namespace cath { namespace test { } }

using namespace cath::rslv;
using namespace cath::test;

using boost::make_optional;
using boost::none;
using std::string;

BOOST_AUTO_TEST_SUITE(crh_filter_spec_test_suite)

BOOST_AUTO_TEST_CASE(applies_hmm_coverages_correctly_for_normal_id) {
	const string id = "normal_id";
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{},                                                                      id ), none                 );
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{}.set_min_hmm_coverage_frac( 0.7 ),                                     id ), make_optional( 0.7 ) );
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{}                                 .set_min_dc_hmm_coverage_frac( 0.8 ), id ), none                 );
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{}.set_min_hmm_coverage_frac( 0.7 ).set_min_dc_hmm_coverage_frac( 0.8 ), id ), make_optional( 0.7 ) );
}

BOOST_AUTO_TEST_CASE(applies_hmm_coverages_correctly_for_dc_id) {
	const string id = "dc_6772b88707fa16c1222d117998650d3d";
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{},                                                                      id ), none                 );
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{}.set_min_hmm_coverage_frac( 0.7 ),                                     id ), make_optional( 0.7 ) );
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{}                                 .set_min_dc_hmm_coverage_frac( 0.8 ), id ), make_optional( 0.8 ) );
	BOOST_CHECK_EQUAL( hmm_coverage_for_match( crh_filter_spec{}.set_min_hmm_coverage_frac( 0.7 ).set_min_dc_hmm_coverage_frac( 0.8 ), id ), make_optional( 0.8 ) );
}

BOOST_AUTO_TEST_SUITE_END()

