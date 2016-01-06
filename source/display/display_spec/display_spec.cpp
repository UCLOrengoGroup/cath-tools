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

#include "display_spec.h"

#include <boost/exception/diagnostic_information.hpp>
#include <boost/log/trivial.hpp>
#include <boost/optional.hpp>

#include "common/c++14/make_unique.h"
#include "display/display_colour/display_colour.h"
#include "display/display_colour/display_colour_list.h"
#include "display/display_colour/display_colour_gradient.h"
#include "display/display_colourer/detail/score_colour_handler.h"
#include "display/display_colourer/display_colourer_alignment.h"
#include "display/display_colourer/display_colourer_consecutive.h"

using namespace cath;
using namespace cath::common;
using namespace cath::detail;
using namespace std;

using boost::none;

/// \brief A string value to use internally to indicate colours haven't been specified
///
/// This is used (rather than, say, making display_colours_string an optional<string>) so that
/// display_colours_string can passed to Boost program_options and handled correctly and so that
/// unspecified can easily be distinguished from specified as empty.
const string display_spec::COLOURS_UNSPECIFIED = "this string is used internally to indicate the colours haven't been specified";

/// \brief Ctor for display_spec
display_spec::display_spec(const string &arg_display_colours_string,    ///< TODOCUMENT
                           const bool   &arg_gradient_colour_alignment, ///< TODOCUMENT
                           const bool   &arg_show_scores_if_present,    ///< TODOCUMENT
                           const bool   &arg_scores_to_equivs,          ///< TODOCUMENT
                           const bool   &arg_normalise_scores           ///< TODOCUMENT
                           ) : display_colours_string    ( arg_display_colours_string    ),
                               gradient_colour_alignment ( arg_gradient_colour_alignment ),
                               show_scores_if_present    ( arg_show_scores_if_present    ),
                               scores_to_equivs          ( arg_scores_to_equivs          ),
                               normalise_scores          ( arg_normalise_scores          ) {
}

/// \brief Private getter of whether a display_colours_string has been set
bool display_spec::has_display_colours_string() const {
	return display_colours_string != COLOURS_UNSPECIFIED;
}

/// \brief Private getter of display_colours_string or default value if none has been specified
const string & display_spec::get_display_colours_string_or_default() const {
	return has_display_colours_string() ? display_colours_string
	                                    : display_colour_list::DEFAULT_COLOURS_STRING;
}

/// \brief TODOCUMENT
display_colour_list display_spec::get_colour_list() const {
	return make_display_colour_list_from_string(
		get_display_colours_string_or_default()
	);
}

/// \brief TODOCUMENT
bool display_spec::get_gradient_colour_alignment() const {
	return gradient_colour_alignment;
}

/// \brief TODOCUMENT
bool display_spec::get_show_scores_if_present() const {
	return show_scores_if_present;
}

/// \brief TODOCUMENT
bool display_spec::get_scores_to_equivs() const {
	return scores_to_equivs;
}

/// \brief TODOCUMENT
bool display_spec::get_normalise_scores() const {
	return normalise_scores;
}

/// \brief TODOCUMENT
void display_spec::set_display_colours_string(const string &arg_display_colours_string ///< TODOCUMENT
                                              ) {
	display_colours_string = arg_display_colours_string;
}

/// \brief TODOCUMENT
void display_spec::set_gradient_colour_alignment(const bool &arg_gradient_colour_alignment ///< TODOCUMENT
                                                 ) {
	gradient_colour_alignment = arg_gradient_colour_alignment;
}

/// \brief TODOCUMENT
void display_spec::set_show_scores_if_present(const bool &arg_show_scores_if_present ///< TODOCUMENT
                                              ) {
	show_scores_if_present = arg_show_scores_if_present;
}

/// \brief TODOCUMENT
void display_spec::set_scores_to_equivs(const bool &arg_scores_to_equivs ///< TODOCUMENT
                                        ) {
	scores_to_equivs = arg_scores_to_equivs;
}

/// \brief TODOCUMENT
void display_spec::set_normalise_scores(const bool &arg_normalise_scores ///< TODOCUMENT
                                        ) {
	normalise_scores = arg_normalise_scores;
}

/// \brief TODOCUMENT
unique_ptr<const display_colourer> cath::get_display_colourer(const display_spec            &arg_display_spec,   ///< TODOCUMENT
                                                              const display_colour_gradient &arg_colour_gradient ///< TODOCUMENT
                                                              ) {
	const score_colour_handler colour_handler{
		arg_display_spec.get_show_scores_if_present(),
		arg_display_spec.get_scores_to_equivs(),
		arg_display_spec.get_normalise_scores()
	};
	if ( arg_display_spec.get_gradient_colour_alignment() ) {
		return { common::make_unique< display_colourer_alignment   >( colour_handler, arg_colour_gradient                ) };
	}
	else {
		return { common::make_unique< display_colourer_consecutive >( colour_handler, arg_display_spec.get_colour_list() ) };
	}
}


/// \brief String describing any problems, or "" if none (as part of the interface to friend display_options_block)
opt_str cath::invalid_string(const display_spec &arg_display_spec ///< TODOCUMENT
                             ) {
	try {
		arg_display_spec.get_colour_list();
	}
	catch (const boost::exception &e) {
		return "Colour list could not be parsed from \"" + arg_display_spec.get_display_colours_string_or_default()
			+ "\". Specific error was: "                 + diagnostic_information(e);
	}
	catch (...) {
		return "Colour list could not be parsed from \"" + arg_display_spec.get_display_colours_string_or_default() + "\"";
	}
	return none;
}
