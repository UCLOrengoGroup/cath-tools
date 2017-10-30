/// \file
/// \brief The display_colour_gradient test suite

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

#include "common/boost_addenda/test/boost_check_no_throw_diag.hpp"
#include "display_colour/display_colour.hpp"
#include "display_colour/display_colour_gradient.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

namespace cath {
	namespace test {

		/// \brief The display_colour_gradient_test_suite_fixture to assist in testing display_colour_gradient
		struct display_colour_gradient_test_suite_fixture {
		protected:
			~display_colour_gradient_test_suite_fixture() noexcept = default;
		};

	}
}  // namespace cath

BOOST_FIXTURE_TEST_SUITE(display_colour_gradient_test_suite, cath::test::display_colour_gradient_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ctor_works) {
	BOOST_CHECK_NO_THROW_DIAG( display_colour_gradient( { display_colour::RED, display_colour::BLUE }, 2 ) );
	BOOST_CHECK_NO_THROW_DIAG( display_colour_gradient( { display_colour::RED, display_colour::BLUE }, 1 ) );
	BOOST_CHECK_NO_THROW_DIAG( display_colour_gradient( { display_colour::RED, display_colour::BLUE }, 0 ) );

	BOOST_CHECK_NO_THROW_DIAG( display_colour_gradient( { display_colour::RED                       }, 2 ) );
	BOOST_CHECK_NO_THROW_DIAG( display_colour_gradient( { display_colour::RED                       }, 1 ) );
	BOOST_CHECK_NO_THROW_DIAG( display_colour_gradient( { display_colour::RED                       }, 0 ) );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(ctor_throws_with_no_colours) {
	BOOST_CHECK_THROW( display_colour_gradient( display_colour_vec(), 2 ), invalid_argument_exception );
	BOOST_CHECK_THROW( display_colour_gradient( display_colour_vec(), 1 ), invalid_argument_exception );
	BOOST_CHECK_THROW( display_colour_gradient( display_colour_vec(), 0 ), invalid_argument_exception );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(default_has_correct_properties) {
	constexpr size_t num_steps   = 63                ;
	constexpr size_t half_way    =  1 + num_steps / 2;
	constexpr size_t num_points  =  5                ;
	constexpr size_t num_colours = ( num_steps + 1 ) * ( num_points - 1 ) + 1;

	const display_colour_gradient default_colour_gradient = make_default_colour_gradient();

	BOOST_CHECK_EQUAL( num_steps,                    default_colour_gradient.get_steps_in_between_points()                            );
	BOOST_CHECK_EQUAL( num_points,                   default_colour_gradient.get_num_colour_points()                                  );
	BOOST_CHECK_EQUAL( num_colours,                  get_num_colours(default_colour_gradient)                                         );

	BOOST_CHECK_THROW( get_colour_of_index(default_colour_gradient, num_colours), invalid_argument_exception                          );

	BOOST_CHECK_EQUAL( display_colour::BLUE,          get_colour_of_index( default_colour_gradient, 0 * ( num_steps + 1 )            ) );
	BOOST_CHECK_EQUAL( display_colour(0.0, 0.5, 1.0), get_colour_of_index( default_colour_gradient, 0 * ( num_steps + 1 ) + half_way ) );
	BOOST_CHECK_EQUAL( display_colour::CYAN,          get_colour_of_index( default_colour_gradient, 1 * ( num_steps + 1 )            ) );
	BOOST_CHECK_EQUAL( display_colour(0.0, 1.0, 0.5), get_colour_of_index( default_colour_gradient, 1 * ( num_steps + 1 ) + half_way ) );
	BOOST_CHECK_EQUAL( display_colour::GREEN,         get_colour_of_index( default_colour_gradient, 2 * ( num_steps + 1 )            ) );
	BOOST_CHECK_EQUAL( display_colour(0.5, 1.0, 0.0), get_colour_of_index( default_colour_gradient, 2 * ( num_steps + 1 ) + half_way ) );
	BOOST_CHECK_EQUAL( display_colour::YELLOW,        get_colour_of_index( default_colour_gradient, 3 * ( num_steps + 1 )            ) );
	BOOST_CHECK_EQUAL( display_colour(1.0, 0.5, 0.0), get_colour_of_index( default_colour_gradient, 3 * ( num_steps + 1 ) + half_way ) );
	BOOST_CHECK_EQUAL( display_colour::RED,           get_colour_of_index( default_colour_gradient, 4 * ( num_steps + 1 )            ) );

//	for (const size_t &colour_ctr : indices( num_colours ) ) {
//		cerr << colour_ctr << "\t" << get_colour_of_index(default_colour_gradient, colour_ctr) << endl;
//	}
}

//display_colour_gradient

BOOST_AUTO_TEST_SUITE_END()


