///// \file
///// \brief The superposition I/O test suite

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

#include "superposition/io/superposition_io.hpp"
#include "test/superposition_fixture.hpp"

using namespace cath::geom;
using namespace cath::sup;
using namespace std;

/// \brief The superposition_io_test_suite_fixture to assist in testing superposition I/O
struct superposition_io_test_suite_fixture : protected superposition_fixture{
protected:
	~superposition_io_test_suite_fixture() noexcept = default;
};

BOOST_FIXTURE_TEST_SUITE(superposition_io_test_suite, superposition_io_test_suite_fixture)

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_SUITE(write)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup) {
	BOOST_CHECK_EQUAL( to_json_string( the_sup, false ), sup_json_str );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(write)

BOOST_AUTO_TEST_CASE(from_json_string_works_for_identity) {
	BOOST_CHECK_EQUAL( superposition_from_json_string( sup_json_str ), the_sup );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

