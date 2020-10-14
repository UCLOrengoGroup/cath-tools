/// \file
/// \brief The display_colour class definitions

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

#include "display_colour.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"

#include <tuple>

using namespace cath;
using namespace cath::common;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::join;
using boost::lexical_cast;
using boost::numeric_cast;
using std::tie;

/// \brief TODOCUMENT
const string        display_colour::COMPONENT_SEPARATOR(",");

/// \brief TODOCUMENT
const display_colour display_colour::BLACK        ( 0.00, 0.00, 0.00 );

/// \brief TODOCUMENT
const display_colour display_colour::WHITE        ( 1.00, 1.00, 1.00 );

/// \brief TODOCUMENT
const display_colour display_colour::RED          ( 1.00, 0.00, 0.00 );

/// \brief TODOCUMENT
const display_colour display_colour::LIGHT_RED    ( 1.00, 0.30, 0.30 );

/// \brief TODOCUMENT
const display_colour display_colour::YELLOW       ( 1.00, 1.00, 0.00 );

/// \brief TODOCUMENT
const display_colour display_colour::DARK_YELLOW  ( 0.70, 0.70, 0.00 );

/// \brief TODOCUMENT
const display_colour display_colour::GREEN        ( 0.00, 1.00, 0.00 );

/// \brief TODOCUMENT
const display_colour display_colour::LIGHT_GREEN  ( 0.30, 1.00, 0.30 );

/// \brief TODOCUMENT
const display_colour display_colour::CYAN         ( 0.00, 1.00, 1.00 );

/// \brief TODOCUMENT
const display_colour display_colour::DARK_CYAN    ( 0.00, 0.70, 0.70 );

/// \brief TODOCUMENT
const display_colour display_colour::BLUE         ( 0.00, 0.00, 1.00 );

/// \brief TODOCUMENT
const display_colour display_colour::LIGHT_BLUE   ( 0.30, 0.30, 1.00 );

/// \brief TODOCUMENT
const display_colour display_colour::MAGENTA      ( 1.00, 0.00, 1.00 );

/// \brief TODOCUMENT
const display_colour display_colour::DARK_MAGENTA ( 0.70, 0.00, 0.70 );

/// \brief TODOCUMENT
void display_colour::check_component_value(const double &prm_component_value ///< TODOCUMENT
                                           ) {
	using boost::math::isfinite;
	if ( ! isfinite( prm_component_value ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Viewer colour component value "
			+ lexical_cast<string>(prm_component_value)
			+ " is not a finite number"
		));
	}
	if ( prm_component_value < 0.0 || prm_component_value > 1.0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Viewer colour component value "
			+ lexical_cast<string>( prm_component_value )
			+ " is not between 0 and 1"
		));
	}
}

/// \brief Ctor for display_colour
display_colour::display_colour(const double &prm_r, ///< TODOCUMENT
                               const double &prm_g, ///< TODOCUMENT
                               const double &prm_b  ///< TODOCUMENT
                               ) : r ( prm_r ),
                                   g ( prm_g ),
                                   b ( prm_b ) {
	check_component_value( r );
	check_component_value( g );
	check_component_value( b );
}

/// \brief Getter for the red component (between 0 and 1)
const double & display_colour::get_r() const {
	return r;
}

/// \brief Getter for the green component (between 0 and 1)
const double & display_colour::get_g() const {
	return g;
}

/// \brief Getter for the blue component (between 0 and 1)
const double & display_colour::get_b() const {
	return b;
}

/// \brief Return a copy of the specified colour that has been lightened by the specified fraction
///
/// \relates display_colour
display_colour cath::lighten_by_fraction(const display_colour &prm_colour,  ///< The colour from which to start
                                         const double         &prm_fraction ///< The fraction change (0.0 gives the original colour; 1.0 gives white)
                                         ) {
	return rgb_mid_point( prm_colour, display_colour::WHITE, prm_fraction );
}

/// \brief Return a copy of the specified colour that has been darkened by the specified fraction
///
/// \relates display_colour
display_colour cath::darken_by_fraction(const display_colour &prm_colour,  ///< The colour from which to start
                                        const double         &prm_fraction ///< The fraction change (0.0 gives the original colour; 1.0 gives black)
                                        ) {
	return rgb_mid_point( prm_colour, display_colour::BLACK, prm_fraction );
}

/// \brief Return whether the first colour is less than the second
///
/// This is a simple ordering that compares on the red, green and then blue components
///
/// \todo Is this actually useful? If so, document the context; if not, consider removing.
///
/// \relates display_colour
bool cath::operator<(const display_colour &prm_lhs, ///< The first  display_colour to compare
                     const display_colour &prm_rhs  ///< The second display_colour to compare
                     ) {
	return (
		tie( prm_lhs.get_r(), prm_lhs.get_g(), prm_lhs.get_b() )
		<
		tie( prm_rhs.get_r(), prm_rhs.get_g(), prm_rhs.get_b() )
	);
}

/// \brief TODOCUMENT
///
/// \relates display_colour
display_colour cath::display_colour_from_string(const string &prm_colour ///< TODOCUMENT
                                                ) {
	const str_vec component_strings = split_build<str_vec>( prm_colour, is_any_of( display_colour::COMPONENT_SEPARATOR ) );
	doub_vec component_numbers;

	component_numbers.reserve( component_strings.size() );
	for (const string &component_string : component_strings) {
		// Using lexical_cast_stringstream rather than lexical_cast because lexical_cast<double>() seems
		// to (sometimes?) have problems when run under valgrind in that it converts "1.000" to a number
		// that's a tiny bit bigger than 1 (visible using std::setprecision(50))
		//
		// Might relate to this: https://bugzilla.redhat.com/show_bug.cgi?id=837650 ?
		//
		// ...now trying to switch to stod instead
		//
		/// \todo If stod() appears to be an adeqate replacement, then drop lexical_cast_stringstream()
		component_numbers.push_back( stod( component_string ) );
	}
	return display_colour_from_components(component_numbers);
}

/// \brief TODOCUMENT
///
/// \relates display_colour
display_colour cath::display_colour_from_components(const doub_vec &prm_components ///< TODOCUMENT
                                                    ) {
	if (prm_components.size() != 3) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct viewer colour from component values because there aren't three of them"));
	}
	return display_colour(
		prm_components[ 0 ],
		prm_components[ 1 ],
		prm_components[ 2 ]
	);
}

/// \brief TODOCUMENT
///
/// \relates display_colour
string cath::hex_string_of_colour(const display_colour &prm_display_colour ///< TODOCUMENT
                                  ) {
	ostringstream output_ss;
	const doub_vec comps = { prm_display_colour.get_r(),
	                         prm_display_colour.get_g(),
	                         prm_display_colour.get_b() };
	for (const double &comp : comps) {
		output_ss << hex << setfill('0') << setw(2) << numeric_cast<size_t>( 255 * comp );
	}
	return output_ss.str();
}

/// \brief TODOCUMENT
///
/// \relates display_colour
string cath::comma_separated_string_of_display_colour(const display_colour &prm_display_colour ///< TODOCUMENT
                                                      ) {
	const str_vec component_strings = {
		lexical_cast<string>( prm_display_colour.get_r() ),
		lexical_cast<string>( prm_display_colour.get_g() ),
		lexical_cast<string>( prm_display_colour.get_b() )
	};
	return join( component_strings, display_colour::COMPONENT_SEPARATOR );
}

/// \brief TODOCUMENT
///
/// \relates display_colour
display_colour cath::rgb_mid_point(const display_colour &prm_colour1,      ///< TODOCUMENT
                                   const display_colour &prm_colour2,      ///< TODOCUMENT
                                   const double         &prm_fraction_thru ///< TODOCUMENT
                                   ) {
	using boost::math::isfinite;
	if ( ! isfinite( prm_fraction_thru ) || prm_fraction_thru < 0.0 || prm_fraction_thru > 1.0 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Fraction through must be a finite number between 0 and 1"));
	}

	const double fraction_from = 1.0 - prm_fraction_thru;
	return display_colour(
		fraction_from * prm_colour1.get_r() + prm_fraction_thru * prm_colour2.get_r(),
		fraction_from * prm_colour1.get_g() + prm_fraction_thru * prm_colour2.get_g(),
		fraction_from * prm_colour1.get_b() + prm_fraction_thru * prm_colour2.get_b()
	);
}

/// \brief TODOCUMENT
///
/// \relates display_colour
ostream & cath::operator<<(ostream              &prm_ostream,      ///< TODOCUMENT
                           const display_colour &prm_display_colour ///< TODOCUMENT
                           ) {
	prm_ostream << "display_colour[" << prm_display_colour.get_r();
	prm_ostream << ", "              << prm_display_colour.get_g();
	prm_ostream << ", "              << prm_display_colour.get_b();
	prm_ostream << "]";
	return prm_ostream;
}

