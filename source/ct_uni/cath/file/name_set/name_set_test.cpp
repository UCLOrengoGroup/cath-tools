/// \file
/// \brief The name_set test suite

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

#include "cath/file/name_set/name_set.hpp"

namespace cath { namespace test { } }

using namespace ::cath::file;
using namespace ::cath::test;
using namespace ::std::literals::string_literals;

namespace cath {
	namespace test {

		/// \brief The name_set_test_suite_fixture to assist in testing name_set
		struct name_set_test_suite_fixture {
		protected:
			~name_set_test_suite_fixture() noexcept = default;
		};

	} // namespace test
} // namespace cath

BOOST_FIXTURE_TEST_SUITE(name_set_test_suite, name_set_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_CHECK_EQUAL( to_string( name_set( "a"s             ) ), R"(name_set[name_from_acq:"a"])"                          );
	BOOST_CHECK_EQUAL( to_string( name_set( "a"s, "b"s       ) ), R"(name_set[name_from_acq:"a", spec_id:"b"])"             );
	BOOST_CHECK_EQUAL( to_string( name_set( "a"s, "b"s, "c"s ) ), R"(name_set[name_from_acq:"a", spec_id:"b", dom_id:"c"])" );
}

BOOST_AUTO_TEST_SUITE_END()

