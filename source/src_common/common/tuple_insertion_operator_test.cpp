/// \file
/// \brief The tuple_insertion_operator test suite

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

#include <boost/lexical_cast.hpp>
#include <boost/test/auto_unit_test.hpp>

#include "tuple_insertion_operator.hpp"

#include "common/size_t_literal.hpp"

using namespace cath::common;

using boost::lexical_cast;
using std::make_tuple;
using std::string;

BOOST_AUTO_TEST_SUITE(tuple_insertion_operator_test_suite)

BOOST_AUTO_TEST_CASE(basic) {
	BOOST_TEST( tuple_to_string     ( make_tuple( 3_z, 2.3 ) ) == "std::tuple<unsigned long, double>(3, 2.300000)" );
	BOOST_TEST( lexical_cast<string>( make_tuple( 3_z, 2.3 ) ) == "std::tuple<unsigned long, double>(3, 2.300000)" );
}

BOOST_AUTO_TEST_SUITE_END()
