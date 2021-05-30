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
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"

#include <tuple>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

using ::boost::algorithm::is_any_of;
using ::boost::algorithm::join;
using ::boost::lexical_cast;
using ::boost::numeric_cast;

/// \brief TODOCUMENT
///
/// \relates display_colour
display_colour cath::display_colour_from_string(const string &prm_colour ///< TODOCUMENT
                                                ) {
	const auto component_strings = split_build<str_vec>( prm_colour, is_any_of( display_colour::COMPONENT_SEPARATOR ) );
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
		/// \todo If stod() appears to be an adequate replacement, then drop lexical_cast_stringstream()
		component_numbers.push_back( stod( component_string ) );
	}

	if ( component_numbers.size() != 3 ) {
		BOOST_THROW_EXCEPTION( invalid_argument_exception(
		  "Cannot construct display_colour from component values from string because there aren't three of them" ) );
	}

	return display_colour_from_components( make_array_of_first_n<3>( component_numbers ) );
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
ostream & cath::operator<<(ostream              &prm_ostream,      ///< TODOCUMENT
                           const display_colour &prm_display_colour ///< TODOCUMENT
                           ) {
	prm_ostream << "display_colour[" << prm_display_colour.get_r();
	prm_ostream << ", "              << prm_display_colour.get_g();
	prm_ostream << ", "              << prm_display_colour.get_b();
	prm_ostream << "]";
	return prm_ostream;
}

