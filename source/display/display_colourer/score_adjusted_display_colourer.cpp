/// \file
/// \brief The score_adjusted_display_colourer class definitions

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

#include "score_adjusted_display_colourer.h"

#include "alignment/alignment_context.h"
#include "display/display_colourer/display_colour_spec.h"
#include "exception/invalid_argument_exception.h"

using namespace cath;
using namespace cath::align;
using namespace cath::detail;
// using namespace std;

/// \brief TODOCUMENT
display_colour_spec score_adjusted_display_colourer::do_get_colour_spec(const alignment_context &arg_alignment_context ///< TODOCUMENT
                                                                        ) const {
	display_colour_spec        unscored_spec = do_get_colour_spec_before_scoring( arg_alignment_context );
	const alignment           &the_alignment = arg_alignment_context.get_alignment();
	const alignment::size_type num_entries   = the_alignment.num_entries();
	const alignment::size_type aln_length    = the_alignment.length();

	if ( get_score_handler().get_show_scores_if_present() ) {
		for (alignment::size_type entry = 0; entry < num_entries; ++entry) {
			for (size_t index = 0; index < aln_length; ++index) {
				const aln_posn_opt position = the_alignment.position_of_entry_of_index( entry, index );
				if ( position ) {
					const display_colour_opt base_col = unscored_spec.get_base_clr();
					const display_colour_opt pdb_col  = get_clr_of_pdb_index          ( unscored_spec, entry            );
					const display_colour_opt res_col  = get_clr_of_pdb_and_res_indices( unscored_spec, entry, *position );
					const display_colour best_colour = res_col  ? ( *res_col  ) :
					                                   pdb_col  ? ( *pdb_col  ) :
					                                   base_col ? ( *base_col ) : display_colour::BLACK;
					const display_colour the_colour = score_colour_copy( get_score_handler(), the_alignment, entry, index, best_colour );
					unscored_spec.colour_pdb_residue(
						entry,
						*position,
						the_colour,
						true
					);
				}
			}
		}
	}

	return unscored_spec;
}

/// \brief TODOCUMENT
const score_colour_handler & score_adjusted_display_colourer::get_score_handler() const {
	return score_handler;
}

/// \brief Ctor for score_adjusted_display_colourer
score_adjusted_display_colourer::score_adjusted_display_colourer(const score_colour_handler &arg_score_handler ///< TODOCUMENT
                                                                 ) : score_handler( arg_score_handler ) {
}
