/// \file
/// \brief The hits_processor_list test suite

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

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/test/unit_test.hpp>

#include "cath/resolve_hits/read_and_process_hits/hits_processor/hits_processor_list.hpp"
#include "cath/resolve_hits/read_and_process_hits/hits_processor/summarise_hits_processor.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

#include <sstream>

namespace cath { namespace test { } }

using namespace ::cath::rslv;
using namespace ::cath::rslv::detail;
using namespace ::cath::test;

using ::boost::algorithm::all_of;
using ::std::make_unique;
using ::std::ostream;
using ::std::ostringstream;

namespace cath {
	namespace test {

		/// \brief The hits_processor_list_test_suite_fixture to assist in testing hits_processor_list
		struct hits_processor_list_test_suite_fixture {
		protected:
			~hits_processor_list_test_suite_fixture() noexcept = default;

			/// \brief A dummy ostringstream to which tests can write
			ostringstream test_ss;

			/// \brief A dummy ostream_ref to which the tests can write
			ostream_ref test_ss_ref{ test_ss };

			/// \brief A dummy ostream_ref_vec to which the tests can write
			ref_vec<ostream> test_ss_ref_vec{ test_ss_ref };

			/// \brief A dummy crh_score_spec to be used in tests
			const crh_score_spec score_spec{};

			/// \brief A dummy crh_segment_spec to be used in tests
			const crh_segment_spec segment_spec{};
		};

	}  // namespace test
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(hits_processor_list_test_suite, hits_processor_list_test_suite_fixture)

BOOST_AUTO_TEST_CASE(ctor_does_not_throw) {
	BOOST_CHECK_NO_THROW_DIAG( hits_processor_list{} );
}

BOOST_AUTO_TEST_CASE(size_and_empty_work) {
	BOOST_CHECK      (   hits_processor_list{}.empty() );
	BOOST_CHECK      ( ! hits_processor_list{}.add_processor( summarise_hits_processor             { test_ss_ref_vec } ).empty() );
	BOOST_CHECK      ( ! hits_processor_list{}.add_processor( make_unique<summarise_hits_processor>( test_ss_ref_vec ) ).empty() );
	BOOST_CHECK_EQUAL(   hits_processor_list{}.size(), 0 );
	BOOST_CHECK_EQUAL(   hits_processor_list{}.add_processor( summarise_hits_processor             { test_ss_ref_vec } ).size(), 1 );
	BOOST_CHECK_EQUAL(   hits_processor_list{}.add_processor( make_unique<summarise_hits_processor>( test_ss_ref_vec ) ).size(), 1 );
}

BOOST_AUTO_TEST_CASE(accesses_work) {
	const bool correct_parse_hits = summarise_hits_processor{ test_ss_ref_vec }.wants_hits_that_fail_score_filter();
	const auto the_list = [&] {
		hits_processor_list temp_list;
		temp_list
			.add_processor( summarise_hits_processor{ test_ss_ref_vec } )
			.add_processor( summarise_hits_processor{ test_ss_ref_vec } );
		return temp_list;
	} ();

	BOOST_CHECK_EQUAL( the_list.size(), 2 );
	BOOST_CHECK_EQUAL( the_list[ 0 ].wants_hits_that_fail_score_filter(), correct_parse_hits );
	BOOST_CHECK_EQUAL( the_list[ 1 ].wants_hits_that_fail_score_filter(), correct_parse_hits );
	
	BOOST_CHECK( all_of( the_list, [&] (const hits_processor &x) { return x.wants_hits_that_fail_score_filter() == correct_parse_hits; } ) );
}

BOOST_AUTO_TEST_SUITE_END()
