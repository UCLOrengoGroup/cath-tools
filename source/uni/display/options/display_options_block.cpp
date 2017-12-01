/// \file
/// \brief The display_options_block class definitions

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

#include "display_options_block.hpp"

#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.hpp"
#include "display_colour/display_colour_list.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;

/// \brief The option name for the string to specify the colours to use
const string display_options_block::PO_VIEWER_COLOURS            { "viewer-colours"            };

/// \brief The option name for whether to display a gradient of colours
const string display_options_block::PO_GRADIENT_COLOUR_ALIGNMENT { "gradient-colour-alignment" };

/// \brief The option name for whether to use colour to indicate scores (if they're present)
const string display_options_block::PO_SHOW_SCORES_IF_PRESENT    { "show-scores-if-present"    };

/// \brief The option name for whether to colour based on scores to the *present* equivalent positions
const string display_options_block::PO_SCORES_TO_EQUIVS          { "scores-to-equivs"          };

/// \brief The option name for whether to colour based on scores normalised across the alignment, rather than absolute scores
const string display_options_block::PO_NORMALISE_SCORES          { "normalise-scores"          };

/// \brief A standard do_clone method
unique_ptr<options_block> display_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string display_options_block::do_get_block_name() const {
	return "Viewer (eg PyMOL, Jmol etc) options";
}

/// \brief Add this block's options to the provided options_description
void display_options_block::do_add_visible_options_to_description(options_description &arg_desc,           ///< The options_description to which the options are added
                                                                  const size_t        &/*arg_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                  ) {
	const string colrs_varname{ "<colrs>" };

	const auto display_colours_string_notifier    = [&] (const string &x) { the_display_spec.set_display_colours_string   ( x ); };
	const auto gradient_colour_alignment_notifier = [&] (const bool   &x) { the_display_spec.set_gradient_colour_alignment( x ); };
	const auto show_scores_if_present_notifier    = [&] (const bool   &x) { the_display_spec.set_show_scores_if_present   ( x ); };
	const auto scores_to_equivs_notifier          = [&] (const bool   &x) { the_display_spec.set_scores_to_equivs         ( x ); };
	const auto normalise_scores_notifier          = [&] (const bool   &x) { the_display_spec.set_normalise_scores         ( x ); };
	arg_desc.add_options()
		(
			PO_VIEWER_COLOURS.c_str(),
			value<string>()
				->value_name( colrs_varname                   )
				->notifier  ( display_colours_string_notifier ),
			( "Use " + colrs_varname + " to colour successive entries in the viewer\n"
			  "(format: colon-separated list of comma-separated triples of RGB values between 0 and 1)\n"
			  "(will wrap-around when it runs out of colours)" ).c_str()
		)
		(
			PO_GRADIENT_COLOUR_ALIGNMENT.c_str(),
			bool_switch()
				->notifier     ( gradient_colour_alignment_notifier )
				->default_value( false                              ),
			"Colour the length of the alignment with a rainbow gradient (blue -> red)"
		)
		(
			PO_SHOW_SCORES_IF_PRESENT.c_str(),
			bool_switch()
				->notifier     ( show_scores_if_present_notifier )
				->default_value( false                           ),
			( "Show the alignment scores\n(use with " + PO_GRADIENT_COLOUR_ALIGNMENT + ")" ).c_str()
		)
		(
			PO_SCORES_TO_EQUIVS.c_str(),
			bool_switch()
				->notifier     ( scores_to_equivs_notifier )
				->default_value( false                     ),
			( "Show the alignment scores to equivalent positions, which increases relative scores where few entries are aligned\n"
				"(use with --" + PO_GRADIENT_COLOUR_ALIGNMENT + " and --" + PO_SHOW_SCORES_IF_PRESENT + ")" ).c_str()
		)
		(
			PO_NORMALISE_SCORES.c_str(),
			bool_switch()
				->notifier     ( normalise_scores_notifier )
				->default_value( false                     ),
			( "When showing scores, normalise them to the highest score in the alignment\n"
				"(use with --" + PO_GRADIENT_COLOUR_ALIGNMENT + " and --" + PO_SHOW_SCORES_IF_PRESENT + ")" ).c_str() 
		);
}

/// \brief TODOCUMENT
str_opt display_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                 ) const {
	return ::cath::invalid_string( get_display_spec() );
}

/// \brief Return all options names for this block
str_vec display_options_block::do_get_all_options_names() const {
	return {
		display_options_block::PO_VIEWER_COLOURS,
		display_options_block::PO_GRADIENT_COLOUR_ALIGNMENT,
		display_options_block::PO_SHOW_SCORES_IF_PRESENT,
		display_options_block::PO_SCORES_TO_EQUIVS,
		display_options_block::PO_NORMALISE_SCORES,
	};
}

/// \brief Public, by-value getter for the_display_spec
display_spec display_options_block::get_display_spec() const {
	return the_display_spec;
}
