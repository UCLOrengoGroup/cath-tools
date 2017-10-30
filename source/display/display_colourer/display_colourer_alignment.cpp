/// \file
/// \brief The display_colourer_alignment class definitions

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

#include "display_colourer_alignment.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/numeric.hpp>

#include "alignment/alignment_context.hpp"
#include "alignment/io/outputter/horiz_align_outputter.hpp"
#include "chopping/region/region.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "display/display_colour_spec/display_colour_spec.hpp"
#include "display/viewer/viewer.hpp"
#include "display_colour/display_colour.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "superposition/superposition_context.hpp"

#include <iostream>
#include <numeric>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::sup;

using boost::accumulate;
using boost::lexical_cast;
using std::string;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<display_colourer> display_colourer_alignment::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Implementation of colouring that threads a colour gradient along the length of the alignment,
///        stretching it in proportion to the alignment's residue scores if they're present
///
/// See the class notes for more details
display_colour_spec display_colourer_alignment::do_get_colour_spec(const alignment_context &arg_alignment_context ///< The alignment_context to gradient-colour
                                                                   ) const {
	// Grab some basic details of the alignment and its score
	// (total score is less half the first and half the last scores because they won't be counted
	const alignment           &the_alignment = arg_alignment_context.get_alignment();
	const alignment::size_type num_entries   = the_alignment.num_entries();
	const alignment::size_type aln_length    = the_alignment.length();
	const float_score_vec      scores        = get_total_score_or_num_positions_by_index( the_alignment );
	const double               total_score   = accumulate( scores, 0.0 ) - ( 0.5 * scores.front() )
	                                                                     - ( 0.5 * scores.back()  );

	// Prepare the temporaries to alter while looping over the alignment
	double scores_so_far = 0.0;
	display_colour_spec new_spec;

	new_spec.colour_base( display_colour::BLACK );

	// Loop along the length of the alignment
	for (const size_t &index : indices( aln_length ) ) {
		// If this isn't the first index in the alignment, then increment scores_so_far
		// by the score encountered when moving from the previous residue to this
		// (half the previous residue's score plus half this residue's score)
		if ( index > 0 ) {
			scores_so_far += 0.5 * ( scores[ index - 1 ] + scores[ index ] );
		}

		// Calculate the colour that represents moving through the gradient as much as scores_so_far is through total_score
		const double fraction_through = scores_so_far / total_score;
		const display_colour the_colour = get_colour_of_fraction( gradient, fraction_through );

		// For each entry, determine whether there is a position in the alignment and store the colour there if so
		for (const alignment::size_type &entry : indices( num_entries ) ) {
			const aln_posn_opt position = the_alignment.position_of_entry_of_index( entry, index );
			if ( position ) {
				new_spec.colour_pdb_residue( entry, *position, the_colour );
			}
		}
	}
	return new_spec;
}

/// \brief Return the label for this display_colourer
string display_colourer_alignment::do_get_label() const {
	return "Colour by alignbow";
}

/// \brief Get the gradient with which this display_colourer_alignment should colour alignments
const display_colour_gradient & display_colourer_alignment::get_gradient() const {
	return gradient;
}

/// \brief Ctor
display_colourer_alignment::display_colourer_alignment(display_colour_gradient arg_gradient ///< The gradient with which this display_colourer_alignment should colour alignments
                                                       ) : gradient { std::move( arg_gradient ) } {
}

/// \brief Ctor
display_colourer_alignment::display_colourer_alignment(display_colour_gradient     arg_gradient,     ///< The gradient with which this display_colourer_alignment should colour alignments
                                                       const score_colour_handler &arg_score_handler ///< Specification for post-modifying the colouring based on scores
                                                       ) : super    { arg_score_handler         },
                                                           gradient { std::move( arg_gradient ) } {
}
