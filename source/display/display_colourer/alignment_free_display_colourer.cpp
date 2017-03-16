/// \file
/// \brief The alignment_free_display_colourer class definitions

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

#include "alignment_free_display_colourer.hpp"

#include "alignment/alignment_context.hpp"
#include "display/display_colour_spec/display_colour_spec.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::detail;

/// \brief Implementation of getting the colour spec from the alignment
///
/// For an alignment_free_display_colourer, this just means calling do_get_colour_spec_from_num_entries()
/// with the number of entries
display_colour_spec alignment_free_display_colourer::do_get_colour_spec(const alignment_context &arg_alignment_context ///< The alignment_context to use to colour this alignment
                                                                        ) const {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return display_colour_spec{
		do_get_colour_spec_from_num_entries(
			arg_alignment_context.get_alignment().num_entries()
		)
	};
}

/// \brief Ctor
alignment_free_display_colourer::alignment_free_display_colourer(const score_colour_handler &arg_score_handler ///< Specification for post-modifying the colouring based on scores
                                                                 ) : super( arg_score_handler ) {
}

/// \brief NVI wrapper to get colouring based on the number of entries (structures)
broad_display_colour_spec alignment_free_display_colourer::get_colour_spec_from_num_entries(const size_t &arg_num_entries ///< The number of entries (structures) to colour
                                                                                            ) const {
	return do_get_colour_spec_from_num_entries( arg_num_entries );
}