/// \file
/// \brief The resolve_hits_html_outputter test suite

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

#include <boost/algorithm/string/predicate.hpp>
#include <boost/test/auto_unit_test.hpp>

// #include "common/file/spew.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "resolve_hits/full_hit_list.hpp"
#include "resolve_hits/html_output/resolve_hits_html_outputter.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"
#include "resolve_hits/trim/trim_spec.hpp"

#include <cmath>
#include <iostream>

namespace cath { namespace test { } }

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::test;
using namespace std::literals::string_literals;

using std::log;
using std::pow;

namespace cath {
	namespace test {

		/// \brief The resolve_hits_html_outputter_test_suite_fixture to assist in testing resolve_hits_html_outputter
		struct resolve_hits_html_outputter_test_suite_fixture {
		protected:
			~resolve_hits_html_outputter_test_suite_fixture() noexcept = default;

			/// \brief Make an example calc_hit_list for testing
			full_hit_list make_eg_full_hit_list() {
				
				return full_hit_list{
					{
						full_hit{ { hit_seg_of_res_idcs(   2,  68 ), hit_seg_of_res_idcs( 168, 332 ), }, "1pkyA02",  4.1e-85, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs( 333,                                  464 ), }, "1e0tA01",  2.7e-60, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2,  69 ), hit_seg_of_res_idcs( 167, 336 ), }, "2e28A01",  4.4e-55, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2, 135 ), hit_seg_of_res_idcs( 162, 328 ), }, "3gr4A02",  2.2e-54, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(  69,                                  167 ), }, "1e0tA03",  6.6e-51, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2,  68 ), hit_seg_of_res_idcs( 168, 329 ), }, "3qv9A02",  2.5e-50, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   1,  69 ), hit_seg_of_res_idcs( 167, 328 ), }, "3t05A01",  2.1e-49, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   3,  71 ), hit_seg_of_res_idcs( 167, 329 ), }, "3gg8A02",  4.8e-49, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2, 151 ), hit_seg_of_res_idcs( 168, 326 ), }, "1a3wA02",  1.8e-48, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2,  68 ), hit_seg_of_res_idcs( 161, 329 ), }, "3hqnA02",  3.1e-48, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2,  68 ), hit_seg_of_res_idcs( 169, 329 ), }, "3khdA02",  3.8e-46, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(   2,  71 ), hit_seg_of_res_idcs( 167, 333 ), }, "4drsA02",  5.6e-42, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs( 104,                                  167 ), }, "1pkyC03",  2.2e-30, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs( 170,                                  328 ), }, "3qtgA01",  2.5e-21, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs(  70,                                  166 ), }, "4drsA03",  3.0e-19, hit_score_type::FULL_EVALUE },
						full_hit{ { hit_seg_of_res_idcs( 351,                                  467 ), }, "2e28A03",  3.4e-15, hit_score_type::FULL_EVALUE },
					}
				};
			}

			const full_hit_list eg_full_hit_list = make_eg_full_hit_list();
		};

	}
}

BOOST_FIXTURE_TEST_SUITE(resolve_hits_html_outputter_test_suite, resolve_hits_html_outputter_test_suite_fixture)

BOOST_AUTO_TEST_CASE(step_for_length_works) {
	BOOST_CHECK_THROW( resolve_hits_html_outputter::step_for_length(       0 ), invalid_argument_exception );
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(       1 ),                          1 ); //  1.0
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(       3 ),                          1 ); //  3.0
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(       4 ),                          1 ); //  4.0
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(       5 ),                          1 ); //  5.0
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(       7 ),                          1 ); //  7.0
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(       8 ),                          2 ); //  4.0
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(      19 ),                          2 ); //  9.5
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(      20 ),                          5 ); //  4.0

	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(     200 ),                         50 ); //  4.000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(     399 ),                         50 ); //  7.980
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(     400 ),                        100 ); //  4.000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(     799 ),                        100 ); //  7.990
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(     800 ),                        200 ); //  4.000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(    1999 ),                        200 ); //  9.995
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(    2000 ),                        500 ); //  4.000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(    3999 ),                        500 ); //  7.998

	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(   20000 ),                       5000 ); //  7.9998
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(   39999 ),                       5000 ); //  9.99980
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(   40000 ),                      10000 ); //  4.00000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(   79999 ),                      10000 ); //  7.99990
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(   80000 ),                      20000 ); //  4.00000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(  199999 ),                      20000 ); //  9.99995
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(  200000 ),                      50000 ); //  4.00000
	BOOST_CHECK_EQUAL( resolve_hits_html_outputter::step_for_length(  399999 ),                      50000 ); //  7.99998
}

BOOST_AUTO_TEST_CASE(basic) {
	const auto html_string = resolve_hits_html_outputter::output_html(
		"my_query_sequence_id"s,
		eg_full_hit_list,
		crh_score_spec{
			crh_score_spec::DEFAULT_APPLY_CATH_RULES,
			crh_score_spec::DEFAULT_LONG_DOMAINS_PREFERENCE,
			1.0
		},
		crh_segment_spec{}
	);
	BOOST_CHECK( boost::algorithm::contains( html_string, R"(<html)" ) );

	// spew( "example.cath-resolve-hits.output.html", html_string );
}

BOOST_AUTO_TEST_SUITE_END()
