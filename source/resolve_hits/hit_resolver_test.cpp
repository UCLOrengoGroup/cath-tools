/// \file
/// \brief The hit_resolver test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "common/test_predicate/istreams_equal.hpp"
#include "resolve_hits/hit_resolver.hpp"
#include "resolve_hits/options/spec/crh_spec.hpp"
#include "resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"
#include "resolve_hits/test/resolve_hits_fixture.hpp"

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::rslv;
using namespace cath::rslv::detail;
using namespace cath::test;

using std::istringstream;
using std::ostringstream;
using std::string;
using std::stringstream;

namespace cath {
	namespace test {

		/// \brief The hit_resolver_test_suite_fixture to assist in testing calc_hit_list
		struct hit_resolver_test_suite_fixture : protected resolve_hits_fixture {
		protected:
			~hit_resolver_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(hit_resolver_test_suite, hit_resolver_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	istringstream test_iss{ example_input_raw };
	stringstream test_oss;
	read_and_process_mgr the_read_and_process_mgr = make_read_and_process_mgr(
		test_oss,
		crh_spec{}
			.set_score_spec( make_neutral_score_spec() )
	);
	read_hit_list_from_istream( the_read_and_process_mgr, test_iss, hit_score_type::CRH_SCORE );
	istringstream expected_iss{ example_output };

	BOOST_CHECK_ISTREAMS_EQUAL( test_oss, "got_resolved_hits", expected_iss, "istream2" );
}

BOOST_AUTO_TEST_SUITE_END()
