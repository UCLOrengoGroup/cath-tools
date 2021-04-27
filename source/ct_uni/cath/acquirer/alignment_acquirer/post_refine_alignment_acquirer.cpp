/// \file
/// \brief The post_refine_alignment_acquirer class definitions

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

#include "post_refine_alignment_acquirer.hpp"

#include <spdlog/spdlog.h>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/gap/gap_penalty.hpp"
#include "cath/alignment/refiner/alignment_refiner.hpp"
#include "cath/alignment/residue_score/residue_scorer.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::align::gap;
using namespace ::cath::align;
using namespace ::cath::file;

using ::std::pair;

/// \brief Get the concrete alignment_acquirer to acquire an alignment and then refine and rescore it
pair<alignment, size_size_pair_vec> post_refine_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context &prm_strucs_context, ///< The structures corresponding to the alignment
                                                                                                       const align_refining &prm_align_refining  ///< How much refining should be done to the alignment
                                                                                                       ) const {
	const auto unrefined_aln_n_spntree = do_get_alignment_and_spanning_tree( prm_strucs_context );

	// If no refinement is required, return the result
	if ( prm_align_refining == align_refining::NO ) {
		return unrefined_aln_n_spntree;
	}

	// Or if only light refinement is required and this has three or more entries, warn and then return the result
	if ( prm_align_refining == align_refining::LIGHT && unrefined_aln_n_spntree.first.num_entries() > 3 ) {
		::spdlog::warn( "" );
		return unrefined_aln_n_spntree;
	}

	// Otherwise, refine...
	const auto backbone_complete_strucs_context = strucs_context_of_backbone_complete_region_limited_subset_pdbs(
		prm_strucs_context
	);
	const protein_list proteins                 = build_protein_list( backbone_complete_strucs_context );
	const alignment    refined_alignment        = alignment_refiner().iterate( unrefined_aln_n_spntree.first, proteins, gap_penalty( 50, 0 ) );
	const alignment    scored_refined_alignment = score_alignment_copy( residue_scorer(), refined_alignment, proteins );

	// Return the result
	return make_pair( scored_refined_alignment, unrefined_aln_n_spntree.second );
}
