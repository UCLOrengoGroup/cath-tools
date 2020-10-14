/// \file
/// \brief The common_residue_score_based_selection_policy class definitions

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

#include "common_residue_score_based_selection_policy.hpp"

#include "cath/alignment/alignment.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::std;

size_vec common_residue_score_based_selection_policy::do_select_common_residues(const alignment                    &prm_alignment, ///< TODOCUMENT
                                                                                const vector<alignment::size_type> &prm_indices,   ///< TODOCUMENT
                                                                                const alignment::size_type         &prm_entry_a,   ///< TODOCUMENT
                                                                                const alignment::size_type         &prm_entry_b    ///< TODOCUMENT
                                                                                ) const {
	// Sanity check the alignment is scored
	if ( ! prm_alignment.is_scored() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot use a common_residue_score_based_selection_policy to select common coords for an unscored alignment"));
	}

	// Grab the scores of the positions that are in common between the entries prm_entry_a and prm_entry_b
	// Also keep track of the original indices of these positions to reconstruct them afterwards
	doub_doub_pair_vec scores;
	scores.reserve( prm_indices.size() );
	for (const size_t &index : prm_indices) {
		if ( num_present_positions_of_index( prm_alignment, index ) > 1 ) {
			scores.push_back( make_pair(
				get_score_of_entry_and_index( prm_alignment, prm_entry_a, index ),
				get_score_of_entry_and_index( prm_alignment, prm_entry_b, index )
			) );
		}
	}

	// Return the results from the concrete class's implementation of the do_select_common_residues_with_scores() method
	return do_select_common_residues_with_scores( scores );
}

