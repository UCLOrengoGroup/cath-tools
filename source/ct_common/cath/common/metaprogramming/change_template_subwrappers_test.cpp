/// \file
/// \brief The change_template_subwrappers test suite

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#include "change_template_subwrappers.hpp"

#include <boost/test/unit_test.hpp>

#include <list>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace ::cath::common;

using ::std::is_same_v;
using ::std::list;
using ::std::tuple;
using ::std::vector;

BOOST_AUTO_TEST_SUITE(change_template_subwrappers_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	static_assert(
		is_same_v<
			change_template_subwrappers_t<
				tuple< tuple<int>, list<double>, tuple<char> >,
				vector
			>,
			tuple< vector<int>, vector<double>, vector<char> >
		>,
		"change_template_subwrappers_t should change tuple< tuple<int>, list<double>, tuple<char> > with vector into tuple< vector<int>, vector<double>, vector<char> >"
	);

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
