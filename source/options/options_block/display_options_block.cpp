/// \file
/// \brief The display_options_block class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "display_options_block.h"

#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.h"
#include "display/display_colour/display_colour_list.h"

using namespace boost::filesystem;
using namespace boost::program_options;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

/// \brief TODOCUMENT
const string display_options_block::PO_VIEWER_COLOURS            ( "viewer-colours"            );
const string display_options_block::PO_GRADIENT_COLOUR_ALIGNMENT ( "gradient-colour-alignment" );
const string display_options_block::PO_SHOW_SCORES_IF_PRESENT    ( "show-scores-if-present"    );
const string display_options_block::PO_SCORES_TO_EQUIVS          ( "scores-to-equivs"          );
const string display_options_block::PO_NORMALISE_SCORES          ( "normalise-scores"          );

const str_vec display_options_block::ALL_BLOCK_POS = { display_options_block::PO_VIEWER_COLOURS,
                                                       display_options_block::PO_GRADIENT_COLOUR_ALIGNMENT,
                                                       display_options_block::PO_SHOW_SCORES_IF_PRESENT,
                                                       display_options_block::PO_SCORES_TO_EQUIVS,
                                                       display_options_block::PO_NORMALISE_SCORES };

/// \brief A standard do_clone method.
unique_ptr<options_block> display_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string display_options_block::do_get_block_name() const {
	return "Viewer (eg PyMOL, Jmol etc) options";
}

/// \brief Add this block's options to the provided options_description
void display_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                 ) {
	string &display_colours_string = get_display_spec_ref().get_display_colours_string_ref();
	bool   &colour_alignment       = get_display_spec_ref().get_gradient_colour_alignment_ref();
	bool   &show_scores_if_present = get_display_spec_ref().get_show_scores_if_present_ref();
	bool   &scores_to_equivs       = get_display_spec_ref().get_scores_to_equivs_ref();
	bool   &normalise_scores       = get_display_spec_ref().get_normalise_scores_ref();

	arg_desc.add_options()
		( PO_VIEWER_COLOURS.c_str(),            value<string>( &display_colours_string ),                       "Use arg to colour successive entries in the viewer\n(format: colon-separated list of comma-separated triples of RGB values between 0 and 1)\n(will wrap-around when it runs out of colours)" )
		( PO_GRADIENT_COLOUR_ALIGNMENT.c_str(), bool_switch( &colour_alignment       )->default_value( false ), "Colour the length of the alignment with a rainbow gradient (blue -> red)" )
		( PO_SHOW_SCORES_IF_PRESENT.c_str(),    bool_switch( &show_scores_if_present )->default_value( false ), ( "Show the alignment scores\n(use with " + PO_GRADIENT_COLOUR_ALIGNMENT + ")" ).c_str() )
		( PO_SCORES_TO_EQUIVS.c_str(),          bool_switch( &scores_to_equivs       )->default_value( false ), ( "Show the alignment scores to equivalent positions, which increases relative scores where few entries are aligned\n(use with --" + PO_GRADIENT_COLOUR_ALIGNMENT + " and --" + PO_SHOW_SCORES_IF_PRESENT + ")" ).c_str() )
		( PO_NORMALISE_SCORES.c_str(),          bool_switch( &normalise_scores       )->default_value( false ), ( "When showing scores, normalise them to the highest score in the alignment\n(use with --" + PO_GRADIENT_COLOUR_ALIGNMENT + " and --" + PO_SHOW_SCORES_IF_PRESENT + ")" ).c_str() );
}

/// \brief TODOCUMENT
opt_str display_options_block::do_invalid_string() const {
	return get_display_spec_ref().invalid_string();
}

/// \brief Private, non-const-reference getter for the_display_spec
display_spec & display_options_block::get_display_spec_ref() {
	return the_display_spec;
}

/// \brief Private, const-reference getter for the_display_spec
const display_spec & display_options_block::get_display_spec_ref() const {
	return the_display_spec;
}

/// \brief Public, by-value getter for the_display_spec
display_spec display_options_block::get_display_spec() const {
	return get_display_spec_ref();
}

/// \brief Return whether or not the specified variables_map has specified any of the options in this block
bool display_options_block::has_specified_options(const variables_map &arg_vm ///< The variables_map to examine
                                                  ) const {
	for (const string &block_po : ALL_BLOCK_POS) {
		if ( arg_vm.count( block_po ) && ! arg_vm[ block_po ].defaulted() ) {
			return true;
		}
	}
	return false;
}