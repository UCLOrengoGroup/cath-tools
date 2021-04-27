/// \file
/// \brief The display_spec class definitions

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

#include "display_spec.hpp"

#include <boost/exception/diagnostic_information.hpp>
#include <boost/optional.hpp>

#include "cath/common/cpp14/make_unique.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_list.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::make_optional;
using ::boost::none;
using ::std::string;

/// \brief A string value to use internally to indicate colours haven't been specified
///
/// This is used (rather than, say, making display_colours_string an optional<string>) so that
/// display_colours_string can passed to Boost program_options and handled correctly and so that
/// unspecified can easily be distinguished from specified as empty.
const string display_spec::COLOURS_UNSPECIFIED = "this string is used internally to indicate the colours haven't been specified";

constexpr bool display_spec::DEFAULT_GRADIENT_COLOUR_ALIGNMENT;
constexpr bool display_spec::DEFAULT_SHOW_SCORES_IF_PRESENT;
constexpr bool display_spec::DEFAULT_SCORES_TO_EQUIVS;
constexpr bool display_spec::DEFAULT_NORMALISE_SCORES;

/// \brief Ctor for display_spec
display_spec::display_spec(string      prm_display_colours_string,    ///< TODOCUMENT
                           const bool &prm_gradient_colour_alignment, ///< Whether to display a gradient of colours
                           const bool &prm_show_scores_if_present,    ///< Whether to use colour to indicate scores (if they're present)
                           const bool &prm_scores_to_equivs,          ///< Whether to colour based on scores to the *present* equivalent positions
                           const bool &prm_normalise_scores           ///< Whether to colour based on scores normalised across the alignment, rather than absolute scores
                           ) : display_colours_string    { std::move( prm_display_colours_string ) },
                               gradient_colour_alignment { prm_gradient_colour_alignment           },
                               show_scores_if_present    { prm_show_scores_if_present              },
                               scores_to_equivs          { prm_scores_to_equivs                    },
                               normalise_scores          { prm_normalise_scores                    } {
}

/// \brief TODOCUMENT
str_opt display_spec::get_display_colours_string() const {
	return ( display_colours_string != COLOURS_UNSPECIFIED ) ? ::boost::make_optional( display_colours_string )
	                                                         : none;
}

/// \brief Getter for whether to display a gradient of colours
bool display_spec::get_gradient_colour_alignment() const {
	return gradient_colour_alignment;
}

/// \brief Getter for whether to use colour to indicate scores (if they're present)
bool display_spec::get_show_scores_if_present() const {
	return show_scores_if_present;
}

/// \brief Getter for whether to colour based on scores to the *present* equivalent positions
bool display_spec::get_scores_to_equivs() const {
	return scores_to_equivs;
}

/// \brief Getter for whether to colour based on scores normalised across the alignment, rather than absolute scores
bool display_spec::get_normalise_scores() const {
	return normalise_scores;
}

/// \brief TODOCUMENT
void display_spec::set_display_colours_string(const string &prm_display_colours_string ///< TODOCUMENT
                                              ) {
	display_colours_string = prm_display_colours_string;
}

/// \brief Setter for whether to display a gradient of colours
void display_spec::set_gradient_colour_alignment(const bool &prm_gradient_colour_alignment ///< Whether to display a gradient of colours
                                                 ) {
	gradient_colour_alignment = prm_gradient_colour_alignment;
}

/// \brief Setter for whether to use colour to indicate scores (if they're present)
void display_spec::set_show_scores_if_present(const bool &prm_show_scores_if_present ///< Whether to use colour to indicate scores (if they're present)
                                              ) {
	show_scores_if_present = prm_show_scores_if_present;
}

/// \brief Setter for whether to colour based on scores to the *present* equivalent positions
void display_spec::set_scores_to_equivs(const bool &prm_scores_to_equivs ///< Whether to colour based on scores to the *present* equivalent positions
                                        ) {
	scores_to_equivs = prm_scores_to_equivs;
}

/// \brief Setter for whether to colour based on scores normalised across the alignment, rather than absolute scores
void display_spec::set_normalise_scores(const bool &prm_normalise_scores ///< Whether to colour based on scores normalised across the alignment, rather than absolute scores
                                        ) {
	normalise_scores = prm_normalise_scores;
}

/// \brief Return whether the specified display_spec implies the need for an alignment
bool cath::requires_alignment(const display_spec &prm_display_spec ///< The display_spec to query
                              ) {
	return prm_display_spec.get_gradient_colour_alignment();
}

/// \brief Return whether the specified display_spec implies a display_colourer_consecutive
bool cath::is_consecutive(const display_spec &prm_display_spec ///< The display_spec to query
                          ) {
	return ! prm_display_spec.get_gradient_colour_alignment();
}

/// \brief Private getter of whether a display_colours_string has been set
bool cath::has_display_colours_string(const display_spec &prm_display_spec ///< TODOCUMENT
                                      ) {
	return static_cast<bool>( prm_display_spec.get_display_colours_string() );
}

/// \brief String describing any problems, or "" if none (as part of the interface to friend display_options_block)
str_opt cath::invalid_string(const display_spec &prm_display_spec ///< TODOCUMENT
                             ) {
	try {
		get_colour_list( prm_display_spec );
	}
	catch (const boost::exception &e) {
		return "Colour list could not be parsed from \"" + prm_display_spec.get_display_colours_string().value_or( "<no-colour-list-specified>" )
		     + "\". Specific error was: "                + diagnostic_information( e );
	}
	catch (...) {
		return "Colour list could not be parsed from \"" + prm_display_spec.get_display_colours_string().value_or( "<no-colour-list-specified>" ) + "\"";
	}
	return none;
}

/// \brief TODOCUMENT
///
/// \relates display_colour_list
display_colour_list cath::get_colour_list(const display_spec &prm_display_spec ///< TODOCUMENT
                                          ) {
	const auto cols_str_opt = prm_display_spec.get_display_colours_string();
	return make_display_colour_list_from_string(
		cols_str_opt.value_or( display_colour_list::DEFAULT_COLOURS_STRING )
	);
}
