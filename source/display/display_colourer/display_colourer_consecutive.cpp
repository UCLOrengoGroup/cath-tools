/// \file
/// \brief The display_colourer_consecutive class definitions

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

#include "display_colourer_consecutive.h"

#include "alignment/alignment_context.h"
#include "common/clone/make_uptr_clone.h"
#include "display/display_colourer/display_colour_spec.h"
#include "display/viewer/viewer.h"
#include "display_colour/display_colour.h"
#include "superposition/superposition_context.h"

#include <algorithm>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::sup;
using namespace std;

/// \brief A standard do_clone method
unique_ptr<display_colourer> display_colourer_consecutive::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
display_colour_spec display_colourer_consecutive::do_get_colour_spec_before_scoring(const alignment_context &arg_alignment_context ///< TODOCUMENT
                                                                                    ) const {
	// A const-reference to the alignment and grab the number of entries in it
	const alignment           &the_alignment = arg_alignment_context.get_alignment();
	const alignment::size_type num_entries   = the_alignment.num_entries();

	// Create a new display_colour_spec and populate it for the entries with colours
	display_colour_spec new_spec;
	for (alignment::size_type entry = 0; entry < num_entries; ++entry) {
		const display_colour full_colour = colour_of_mod_index( colours, entry );
		new_spec.colour_pdb( entry, full_colour );
	}

	// Return the generated display_colour_spec
	return new_spec;
}

/// \brief Getter for the list of colours with which the structures should be coloured
const display_colour_list & display_colourer_consecutive::get_colours() const {
	return colours;
}

/// \brief Ctor for display_colourer_consecutive
display_colourer_consecutive::display_colourer_consecutive(const score_colour_handler &arg_score_handler, ///< TODOCUMENT
                                                           const display_colour_list  &arg_colours        ///< The list of colours with which the structures should be coloured
                                                           ) : super   ( arg_score_handler ),
                                                               colours ( arg_colours       ) {
}

