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

#include "score_colour_handler.h"

#include <iostream>
#include <numeric>

#include "alignment/alignment.h"
#include "common/type_aliases.h"
#include "display/display_colour/display_colour.h"

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
