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

#include "cath/alignment/alignment.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/display/display_colour_spec/display_colour_spec.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::detail;
using namespace cath::align;
using namespace std;

/// \brief Ctor for score_colour_handler
score_colour_handler::score_colour_handler(const bool &prm_show_scores_if_present, ///< TODOCUMENT
                                           const bool &prm_scores_to_equivs,       ///< TODOCUMENT
                                           const bool &prm_normalise_scores        ///< TODOCUMENT
                                           ) : show_scores_if_present ( prm_show_scores_if_present ),
                                               scores_to_equivs       ( prm_scores_to_equivs       ),
                                               normalise_scores       ( prm_normalise_scores       ) {
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
float_score_type score_colour_handler::get_score_of_position(const alignment &prm_alignment, ///< TODOCUMENT
                                                             const size_t    &prm_entry,     ///< TODOCUMENT
                                                             const size_t    &prm_index      ///< TODOCUMENT
                                                             ) const {
	const bool using_scores     = show_scores_if_present && prm_alignment.is_scored();
	const bool using_this_score = using_scores && has_score( prm_alignment.get_alignment_residue_scores(), prm_entry, prm_index );
	return using_this_score ? get_score( prm_alignment.get_alignment_residue_scores(), prm_entry, prm_index, ! scores_to_equivs, normalise_scores )
	                        : 1.0;
}

/// \brief TODOCUMENT
///
/// \relates score_colour_handler
void cath::detail::score_colour(const score_colour_handler &prm_score_colour_handler, ///< TODOCUMENT
                                const alignment            &prm_alignment,            ///< TODOCUMENT
                                const size_t               &prm_entry,                ///< TODOCUMENT
                                const size_t               &prm_index,                ///< TODOCUMENT
                                display_colour             &prm_colour                ///< TODOCUMENT
                                ) {
	const float_score_type score = prm_score_colour_handler.get_score_of_position( prm_alignment, prm_entry, prm_index );
	prm_colour = rgb_mid_point( display_colour::WHITE, prm_colour, score );
}

/// \brief TODOCUMENT
///
/// \relates score_colour_handler
display_colour cath::detail::score_colour_copy(const score_colour_handler &prm_score_colour_handler, ///< TODOCUMENT
                                               const alignment            &prm_alignment,            ///< TODOCUMENT
                                               const size_t               &prm_entry,                ///< TODOCUMENT
                                               const size_t               &prm_index,                ///< TODOCUMENT
                                               display_colour              prm_colour                ///< TODOCUMENT
                                               ) {
	score_colour(
		prm_score_colour_handler,
		prm_alignment,
		prm_entry,
		prm_index,
		prm_colour
	);
	return prm_colour;
}

/// \brief Adjust an existing display_colour_spec,
///        based on the specified score_colour_handler and alignment
///
/// \relates score_colour_handler
void cath::detail::adjust_display_colour_spec(display_colour_spec        &prm_colour_spec,          ///< The display_colour_spec to alter
                                              const score_colour_handler &prm_score_colour_handler, ///< The specification for how to adjust the colours
                                              const alignment            &prm_alignment             ///< The alignment to use to adjust the colours
                                              ) {
	const alignment::size_type num_entries   = prm_alignment.num_entries();
	const alignment::size_type aln_length    = prm_alignment.length();

	if ( prm_score_colour_handler.get_show_scores_if_present() ) {
		for (const alignment::size_type &entry : indices( num_entries ) ) {
			for (const size_t &index : indices( aln_length ) ) {
				const aln_posn_opt position = prm_alignment.position_of_entry_of_index( entry, index );
				if ( position ) {
					const display_colour_opt &base_col = get_base_clr( prm_colour_spec );
					const display_colour_opt  pdb_col  = get_clr_of_pdb_index          ( prm_colour_spec, entry            );
					const display_colour_opt  res_col  = get_clr_of_pdb_and_res_indices( prm_colour_spec, entry, *position );
					const display_colour best_colour = res_col  ? ( *res_col  ) :
					                                   pdb_col  ? ( *pdb_col  ) :
					                                   base_col ? ( *base_col ) : display_colour::BLACK;
					const display_colour the_colour = score_colour_copy( prm_score_colour_handler, prm_alignment, entry, index, best_colour );
					prm_colour_spec.colour_pdb_residue(
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
display_colour_spec cath::detail::adjust_display_colour_spec_copy(display_colour_spec         prm_colour_spec,          ///< The display_colour_spec from which a copy should be taken and altered
                                                                  const score_colour_handler &prm_score_colour_handler, ///< The specification for how to adjust the colours
                                                                  const alignment            &prm_alignment             ///< The alignment to use to adjust the colours
                                                                  ) {
	adjust_display_colour_spec( prm_colour_spec, prm_score_colour_handler, prm_alignment );
	return prm_colour_spec;
}
