/// \file
/// \brief The display_colour_list class definitions

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

#include "display_colour_list.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/log/trivial.hpp>
#include <boost/optional.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::std;

using ::boost::algorithm::is_any_of;
using ::boost::algorithm::join;

// This can be refreshed using:
//  awk < /usr/local/svn/source/update/trunk/ddmake/colourlist.txt '{print "\t\"" $2 "," $3 "," $4 "\", // " $1}'
// followed by a bit of tidying up
str_vec display_colour_list::DEFAULT_COLOURS_STRING_PARTS= {
	"0.000,0.000,1.000", // Blue
	"1.000,0.000,0.000", // Red
	"0.000,1.000,0.000", // Green
	"1.000,1.000,0.000", // Yellow
	"1.000,0.396,0.459", // Pink
	"0.500,0.500,0.500", // Grey
	"0.627,0.125,0.941", // Purple
	"0.686,0.839,1.000", // BlueTint
	"0.549,0.941,0.549", // LightGreen
	"1.000,0.647,0.000", // Orange
	"0.000,1.000,1.000", // Cyan
	"0.686,0.459,0.349", // Brown
	"0.180,0.545,0.341", // GreenBlue
	"1.000,0.000,0.396", // HotPink
	"1.000,0.000,1.000", // Magenta
	"1.000,0.671,0.733", // PinkTint
	"0.965,0.965,0.459", // YellowTint
	"1.000,0.612,0.000", // Gold
	"0.596,1.000,0.702", // GreenTint
	"1.000,0.271,0.000", // RedOrange
	"0.000,0.980,0.427", // SeaGreen
	"0.227,0.565,1.000", // SkyBlue
	"0.933,0.510,0.933"  // Violet
};
string display_colour_list::COLOURS_SEPARATOR(":");
string display_colour_list::DEFAULT_COLOURS_STRING(join(DEFAULT_COLOURS_STRING_PARTS, COLOURS_SEPARATOR));

/// \brief Ctor for display_colour_list
display_colour_list::display_colour_list(display_colour_vec prm_colours ///< TODOCUMENT
                                         ) : colours { std::move( prm_colours ) } {
}

/// \brief TODOCUMENT
size_t display_colour_list::size() const {
	return colours.size();
}

/// \brief TODOCUMENT
const display_colour & display_colour_list::colour_of_index(const size_t &prm_index ///< TODOCUMENT
                                                            ) const {
	return colours[prm_index];
}

/// \brief Standard const begin() operator to provide range access
auto display_colour_list::begin() const -> const_iterator {
	return common::cbegin( colours );
}

/// \brief Standard const end() operator to provide range access
auto display_colour_list::end() const -> const_iterator {
	return common::cend( colours );
}

/// \brief TODOCUMENT
///
/// \relates display_colour_list
const display_colour & cath::colour_of_mod_index(const display_colour_list &prm_display_colour_list, ///< TODOCUMENT
                                                 const size_t              &prm_index                ///< TODOCUMENT
                                                 ) {
	const size_t num_colours = prm_display_colour_list.size();
	return prm_display_colour_list.colour_of_index( prm_index % num_colours );
}

/// \brief TODOCUMENT
///
/// \relates display_colour_list
display_colour_list cath::make_display_colour_list_from_string(const string &prm_colours_string ///< TODOCUMENT
                                                               ) {
	const str_vec colour_strings = split_build<str_vec>( prm_colours_string, is_any_of( display_colour_list::COLOURS_SEPARATOR ) );
	display_colour_list result( transform_build<display_colour_vec>(
		colour_strings,
		[] (const string &x) { return display_colour_from_string( x ); }
	) );
	return result;
}

/// \brief Return a display_colour_list of the default list of colours
display_colour_list cath::default_display_colour_list() {
	return make_display_colour_list_from_string( display_colour_list::DEFAULT_COLOURS_STRING );
}