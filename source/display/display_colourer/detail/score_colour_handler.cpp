/// \file
/// \brief The score_colour_handler class definitions

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

#include "score_colour_handler.hpp"

#include <iostream>
#include <numeric>

#include "alignment/alignment.hpp"
#include "common/type_aliases.hpp"
#include "display/display_colour_spec/display_colour_spec.hpp"
#include "display_colour/display_colour.hpp"
#include "display_colour/display_colour_type_aliases.hpp"

using namespace cath;
using namespace cath::detail;
using namespace cath::align;
using namespace std;

/// \brief Ctor for score_colour_handler
score_colour_handler::score_colour_handler(const bool &arg_show_scores_if_present, ///< TODOCUMENT
                                           const bool &arg_scores_to_equivs,       ///< TODOCUMENT
                                           const bool &arg_normalise_scores        ///< TODOCUMENT
                                           ) : show_scores_if_present ( arg_show_scores_if_present ),
                                               scores_to_equivs       ( arg_scores_to_equivs       ),
                                               normalise_scores       ( arg_normalise_scores       ) {
}

/// \brief TODOCUMENT
bool score_colour_handler::get_show_scores_if_present() const {
	return show_scores_if_present;
}

/// \brief TODOCUMENT
bool score_colour_handler::get_scores_to_equivs() const {
	return scores_to_equivs;
}

/// \brief TODOCUMENT
bool score_colour_handler::get_normalise_scores() const {
	return normalise_scores;
}

/// \brief TODOCUMENT
float_score_type score_colour_handler::get_score_of_postion(const alignment &arg_alignment, ///< TODOCUMENT
                                                            const size_t    &arg_entry,     ///< TODOCUMENT
                                                            const size_t    &arg_index      ///< TODOCUMENT
                                                            ) const {
	const bool using_scores     = show_scores_if_present && arg_alignment.is_scored();
	const bool using_this_score = using_scores && has_score( arg_alignment.get_alignment_residue_scores(), arg_entry, arg_index );
	return using_this_score ? get_score( arg_alignment.get_alignment_residue_scores(), arg_entry, arg_index, ! scores_to_equivs, normalise_scores )
	                        : 1.0;
}

/// \brief TODOCUMENT
///
/// \relates score_colour_handler
void cath::detail::score_colour(const score_colour_handler &arg_score_colour_handler, ///< TODOCUMENT
                                const alignment            &arg_alignment,            ///< TODOCUMENT
                                const size_t               &arg_entry,                ///< TODOCUMENT
                                const size_t               &arg_index,                ///< TODOCUMENT
                                display_colour             &arg_colour                ///< TODOCUMENT
                                ) {
	const float_score_type score = arg_score_colour_handler.get_score_of_postion( arg_alignment, arg_entry, arg_index );
	arg_colour = rgb_mid_point( display_colour::WHITE, arg_colour, score );
}

/// \brief TODOCUMENT
///
/// \relates score_colour_handler
display_colour cath::detail::score_colour_copy(const score_colour_handler &arg_score_colour_handler, ///< TODOCUMENT
                                               const alignment            &arg_alignment,            ///< TODOCUMENT
                                               const size_t               &arg_entry,                ///< TODOCUMENT
                                               const size_t               &arg_index,                ///< TODOCUMENT
                                               display_colour              arg_colour                ///< TODOCUMENT
                                               ) {
	score_colour(
		arg_score_colour_handler,
		arg_alignment,
		arg_entry,
		arg_index,
		arg_colour
	);
	return arg_colour;
}

/// \brief Adjust an existing display_colour_spec,
///        based on the specified score_colour_handler and alignment
///
/// \relates score_colour_handler
void cath::detail::adjust_display_colour_spec(display_colour_spec        &arg_colour_spec,          ///< The display_colour_spec to alter
                                              const score_colour_handler &arg_score_colour_handler, ///< The specification for how to adjust the colours
                                              const alignment            &arg_alignment             ///< The alignment to use to adjust the colours
                                              ) {
	const alignment::size_type num_entries   = arg_alignment.num_entries();
	const alignment::size_type aln_length    = arg_alignment.length();

	if ( arg_score_colour_handler.get_show_scores_if_present() ) {
		for (alignment::size_type entry = 0; entry < num_entries; ++entry) {
			for (size_t index = 0; index < aln_length; ++index) {
				const aln_posn_opt position = arg_alignment.position_of_entry_of_index( entry, index );
				if ( position ) {
					const display_colour_opt base_col = get_base_clr( arg_colour_spec );
					const display_colour_opt pdb_col  = get_clr_of_pdb_index          ( arg_colour_spec, entry            );
					const display_colour_opt res_col  = get_clr_of_pdb_and_res_indices( arg_colour_spec, entry, *position );
					const display_colour best_colour = res_col  ? ( *res_col  ) :
					                                   pdb_col  ? ( *pdb_col  ) :
					                                   base_col ? ( *base_col ) : display_colour::BLACK;
					const display_colour the_colour = score_colour_copy( arg_score_colour_handler, arg_alignment, entry, index, best_colour );
					arg_colour_spec.colour_pdb_residue(
						entry,
						*position,
						the_colour,
						true
					);
				}
			}
		}
	}
}

/// \brief Adjust and return a copy of an existing display_colour_spec,
///        based on the specified score_colour_handler and alignment
///
/// \relates score_colour_handler
display_colour_spec cath::detail::adjust_display_colour_spec_copy(display_colour_spec         arg_colour_spec,          ///< The display_colour_spec from which a copy should be taken and altered
                                                                  const score_colour_handler &arg_score_colour_handler, ///< The specification for how to adjust the colours
                                                                  const alignment            &arg_alignment             ///< The alignment to use to adjust the colours
                                                                  ) {
	adjust_display_colour_spec( arg_colour_spec, arg_score_colour_handler, arg_alignment );
	return arg_colour_spec;
}
