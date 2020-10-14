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

#include "display_colourer_consecutive.hpp"


#include "cath/alignment/alignment_context.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/display/display_colour_spec/display_colour_spec.hpp"
#include "cath/display/viewer/viewer.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/superposition/superposition_context.hpp"

#include <algorithm>

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::detail;
using namespace ::cath::sup;

using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<display_colourer> display_colourer_consecutive::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
broad_display_colour_spec display_colourer_consecutive::do_get_colour_spec_from_regions(const region_vec_opt_vec &prm_regions ///< The key regions of the structures
                                                                                        ) const {
	// Create a new display_colour_spec and populate it for the entries with colours
	broad_display_colour_spec new_spec;
	new_spec.colour_base( display_colour::BLACK );
	for (const size_t entry_ctr : indices( prm_regions.size() ) ) {
		// const auto entry_regions = prm_regions[ entry_ctr ];
		// if ( entry_regions ) {
		new_spec.colour_pdb(
			entry_ctr,
			colour_of_mod_index( colours, entry_ctr )
		);
	}

	// Return the generated display_colour_spec
	return new_spec;
}

/// \brief Return the label for this display_colourer
string display_colourer_consecutive::do_get_label() const {
	return "Colour by structure";
}

/// \brief Getter for the list of colours with which the structures should be coloured
const display_colour_list & display_colourer_consecutive::get_colours() const {
	return colours;
}

/// \brief Ctor
display_colourer_consecutive::display_colourer_consecutive(display_colour_list prm_colours ///< The list of colours with which the structures should be coloured
                                                           ) : colours { std::move( prm_colours ) } {
}

/// \brief Ctor
display_colourer_consecutive::display_colourer_consecutive(display_colour_list         prm_colours,      ///< The list of colours with which the structures should be coloured
                                                           const score_colour_handler &prm_score_handler ///< Specification for post-modifying the colouring based on scores
                                                           ) : super   { prm_score_handler        },
                                                               colours { std::move( prm_colours ) } {
}

