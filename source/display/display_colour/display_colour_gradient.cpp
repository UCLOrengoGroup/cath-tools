/// \file
/// \brief The display_colour_gradient class definitions

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

#include "display_colour_gradient.h"

#include <boost/algorithm/clamp.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "exception/invalid_argument_exception.h"
#include "display/display_colour/display_colour.h"

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::numeric_cast;

/// \brief TODOCUMENT
void display_colour_gradient::check_values() const {
	if (colour_points.empty()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Gradient must contain at least one colour"));
	}
}

/// \brief Ctor for display_colour_gradient
display_colour_gradient::display_colour_gradient(const display_colour_vec &arg_colour_points, ///< TODOCUMENT
                                                 const size_t                 &arg_steps_in_between_points ///< TODOCUMENT
                                                 ) : colour_points(arg_colour_points),
                                                     steps_in_between_points(arg_steps_in_between_points) {
	check_values();
}

/// \brief TODOCUMENT
size_t display_colour_gradient::get_steps_in_between_points() const {
	return steps_in_between_points;
}

/// \brief TODOCUMENT
size_t display_colour_gradient::get_num_colour_points() const {
	return colour_points.size();
}

/// \brief TODOCUMENT
const display_colour & display_colour_gradient::get_colour_point_of_index(const size_t &arg_index ///< TODOCUMENT
                                                                          ) const {
	return colour_points[arg_index];
}

/// \brief TODOCUMENT
///
/// \relates display_colour_gradient
size_t cath::get_num_colours(const display_colour_gradient &arg_display_colour_gradient ///< TODOCUMENT
                             ) {
	const size_t steps_in_between_points = arg_display_colour_gradient.get_steps_in_between_points();
	const size_t num_colour_points       = arg_display_colour_gradient.get_num_colour_points();
	return (steps_in_between_points + 1) * (num_colour_points - 1)  + 1;
}

/// \brief TODOCUMENT
///
/// \relates display_colour_gradient
display_colour cath::get_colour_of_fraction(const display_colour_gradient &arg_display_colour_gradient, ///< TODOCUMENT
                                            const double                  &arg_fraction_through         ///< TODOCUMENT
                                            ) {
	if ( arg_fraction_through < 0.0 || arg_fraction_through > 1.0000000001 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot find colour in gradient for fraction that is out of range"));
	}
	const double clamped_fraction_through = clamp( arg_fraction_through, 0.0, 1.0 );
	const size_t num_colours = get_num_colours(arg_display_colour_gradient);
	const size_t colour_index = min(
		num_colours - 1,
		numeric_cast<size_t>( floor( clamped_fraction_through * numeric_cast<double>( num_colours) ) )
	);
	return get_colour_of_index( arg_display_colour_gradient, colour_index );
}

/// \brief TODOCUMENT
///
/// \relates display_colour_gradient
display_colour cath::get_colour_of_index(const display_colour_gradient &arg_display_colour_gradient, ///< TODOCUMENT
                                         const size_t                  &arg_index                    ///< TODOCUMENT
                                         ) {
	const size_t num_colours = get_num_colours(arg_display_colour_gradient);
	if (arg_index >= num_colours) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot get colour of index " + lexical_cast<string>(arg_index)
			+ " because there are only "  + lexical_cast<string>(num_colours)
			+ " colours"
		));
	}

	const size_t steps_in_between_points = arg_display_colour_gradient.get_steps_in_between_points();
	const size_t num_colour_points       = arg_display_colour_gradient.get_num_colour_points();

	const size_t prev_point_index        = arg_index / (steps_in_between_points + 1);
	const size_t steps_after_prev_point  = arg_index % (steps_in_between_points + 1);

	const display_colour &prev_colour   = arg_display_colour_gradient.get_colour_point_of_index(                             prev_point_index       );
	const display_colour &next_colour   = arg_display_colour_gradient.get_colour_point_of_index( min( num_colour_points - 1, prev_point_index + 1 ) );
	const double         fraction_thru = numeric_cast<double>(steps_after_prev_point) / numeric_cast<double>(steps_in_between_points + 1);
	return rgb_mid_point(prev_colour, next_colour, fraction_thru);
}

/// \brief TODOCUMENT
///
/// \relates display_colour_gradient
display_colour_gradient cath::make_default_light_colour_gradient() {
	return display_colour_gradient(
		{
			display_colour::LIGHT_BLUE,
			display_colour::CYAN,
			display_colour::LIGHT_GREEN,
			display_colour::YELLOW,
			display_colour::LIGHT_RED
		},
		63
	);
}

/// \brief TODOCUMENT
///
/// \relates display_colour_gradient
display_colour_gradient cath::make_default_colour_gradient() {
	return display_colour_gradient(
		{
			display_colour::BLUE,
			display_colour::CYAN,
			display_colour::GREEN,
			display_colour::YELLOW,
			display_colour::RED
		},
		63
	);
}

/// \brief TODOCUMENT
///
/// \relates display_colour_gradient
display_colour_gradient cath::make_default_dark_colour_gradient() {
	return display_colour_gradient(
		{
			display_colour::BLUE,
			display_colour::DARK_CYAN,
			display_colour::GREEN,
			display_colour::DARK_YELLOW,
			display_colour::RED
		},
		63
	);
}
