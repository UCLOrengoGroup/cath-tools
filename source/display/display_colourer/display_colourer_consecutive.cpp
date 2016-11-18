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

#include <boost/range/irange.hpp>

#include "alignment/alignment_context.h"
#include "common/clone/make_uptr_clone.h"
#include "common/size_t_literal.h"
#include "display/display_colour_spec/display_colour_spec.h"
#include "display/viewer/viewer.h"
#include "display_colour/display_colour.h"
#include "superposition/superposition_context.h"

#include <algorithm>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::sup;

using boost::irange;
using std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<display_colourer> display_colourer_consecutive::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
broad_display_colour_spec display_colourer_consecutive::do_get_colour_spec_from_num_entries(const size_t &arg_num_entries ///< The number of structures to be coloured
                                                                                            ) const {
	// Create a new display_colour_spec and populate it for the entries with colours
	broad_display_colour_spec new_spec;
	for (const size_t entry_ctr : irange( 0_z, arg_num_entries ) ) {
		new_spec.colour_pdb(
			entry_ctr,
			colour_of_mod_index( colours, entry_ctr )
		);
	}

	// Return the generated display_colour_spec
	return new_spec;
}

/// \brief Getter for the list of colours with which the structures should be coloured
const display_colour_list & display_colourer_consecutive::get_colours() const {
	return colours;
}

/// \brief Ctor
display_colourer_consecutive::display_colourer_consecutive(const display_colour_list &arg_colours ///< The list of colours with which the structures should be coloured
                                                           ) : colours { arg_colours } {
}

/// \brief Ctor
display_colourer_consecutive::display_colourer_consecutive(const display_colour_list  &arg_colours,      ///< The list of colours with which the structures should be coloured
                                                           const score_colour_handler &arg_score_handler ///< Specification for post-modifying the colouring based on scores
                                                           ) : super   { arg_score_handler },
                                                               colours { arg_colours       } {
}

