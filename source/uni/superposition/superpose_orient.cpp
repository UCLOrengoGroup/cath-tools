/// \file
/// \brief The superpose_orient class definitions

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
/// MERCHANTABILITY or ORIENTNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "superpose_orient.hpp"

#include <boost/log/trivial.hpp>

#include "alignment/alignment.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "file/pdb/backbone_complete_indices.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_list.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/geometry/orient.hpp"
#include "superposition/superposition.hpp"

using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;
using namespace cath::sup::detail;

/// \brief For the positions in the specified alignment that pass the specified filter, get their corresponding coordinates in
///        the specified PDBs after transformation by the corresponding part of the specified superposition
template <typename Fn>
inline coord_list cath::sup::detail::get_superposed_filtered_coords_in_aln_order(const superposition &arg_superposition, ///< The superposition by which to transform the coordinates before returning them
                                                                                 const alignment     &arg_alignment,     ///< The alignment from which the positions should be taken
                                                                                 const pdb_list      &arg_pdbs,          ///< The PDBs from which the coordinates should be extracted
                                                                                 Fn                 &&arg_function       ///< The boolean predicate to specify whether a position from the alignment should be included (given the entry and index in the alignment)
                                                                                 ) {
	// Check that the number of entries is the same in the PDBs and alignment
	if ( arg_alignment.num_entries() != arg_pdbs.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot get_superposed_filtered_coords_in_aln_order() for mismatching numbers of entries in the alignment and PDBs"));
	}

	// Get the indices of the backbone_complete residues to use throughout the main loop
	const auto backbone_complete_indices = get_backbone_complete_indices( arg_pdbs );
	for (const size_t &entry : indices( arg_alignment.num_entries() ) ) {
		if ( backbone_complete_indices[ entry ].size() != 1_z + *get_last_present_position_of_entry( arg_alignment, entry ) ) {
			BOOST_LOG_TRIVIAL( warning )
				<< "Whilst getting alignment-ordered coords from alignment/PDBs,"
				<< " found that the number of backbone complete indices in structure "
				<< entry
				<< " is "
				<< backbone_complete_indices[ entry ].size()
				<< " yet the last present position in the alignment in that entry is "
				<< *get_last_present_position_of_entry( arg_alignment, entry );
		}
	}

	coord_list results;
	for (const size_t &index : indices( arg_alignment.length() ) ) {
		for (const size_t &entry : indices( arg_alignment.num_entries() ) ) {
			const auto posn_opt = arg_alignment.position_of_entry_of_index( entry, index );
			if ( posn_opt ) {
				if ( invoke( arg_function, entry, index ) ) {
					results.push_back(
						transform_copy(
							arg_superposition,
							entry,
							get_residue_ca_coord_of_backbone_complete_index(
								arg_pdbs[ entry ],
								backbone_complete_indices[ entry ],
								*posn_opt
							)
						)
					);
				}
			}
		}
	}
	return results;
}

/// \brief Get the list of coordinates to orient for the specified alignment and PDBs,
///        after they're modified by the specified superposition
///
/// The coords are returned in the order in which they appear in the alignment
/// (which is relevant because this ordering is used to select whether to flip
///  the final result on the x, y, z axis or not at all).
coord_list cath::sup::get_coords_to_orient(const superposition &arg_superposition, ///< The superposition with which the structures are to be modified
                                           const alignment     &arg_alignment,     ///< The alignment from which the positions should be extracted
                                           const pdb_list      &arg_pdbs           ///< The PDBs from which the coords should be extracted (and then put through the specified superposition)
                                           ) {
	if ( arg_alignment.num_entries() > 1 ) {
		// First try getting the coords for positions that have alignment scores greater than 0
		if ( arg_alignment.is_scored() ) {
			const alignment_residue_scores &scores = arg_alignment.get_alignment_residue_scores();
			const coord_list results = get_superposed_filtered_coords_in_aln_order(
				arg_superposition,
				arg_alignment,
				arg_pdbs,
				[&] (const size_t &entry, const size_t &index) {
					return get_normalised_score_to_all_entries( scores, entry, index ) > 0;
				}
			);
			if ( ! results.empty() ) {
				return results;
			}
		}

		// Otherwise try getting the coords for positions that have more than one entry aligned
		const coord_list results = get_superposed_filtered_coords_in_aln_order(
			arg_superposition,
			arg_alignment,
			arg_pdbs,
			[&] (const size_t &, const size_t &index) {
				return num_present_positions_of_index( arg_alignment, index ) > 1;
			}
		);
		if ( ! results.empty() ) {
			return results;
		}
	}

	// Otherwise return all coords in alignment order
	return get_superposed_filtered_coords_in_aln_order(
		arg_superposition,
		arg_alignment,
		arg_pdbs,
		[] (const size_t &, const size_t &) { return true; }
	);
}

/// \brief Get the orienting transformation for the specified alignment/pdb_list under the specified superposition
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
coord_rot_pair cath::sup::get_orienting_transformation(const superposition &arg_superposition, ///< The superposition of the structures, that the resulting orientation should follow
                                                       const alignment     &arg_alignment,     ///< The alignment defining the coordinates (including their order) to be oriented
                                                       const pdb_list      &arg_pdbs           ///< The PDBs whose coordinates (after superposing) are to be oriented
                                                       ) {
	return get_orienting_transformation( get_coords_to_orient( arg_superposition, arg_alignment, arg_pdbs ) );
}

/// \brief Modify the specified superposition so as to orient the specified (pre-superposed) coordinates
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
///
/// Note: the order of the coordinates matters
void cath::sup::orient_superposition(superposition    &arg_superposition, ///< The superposition of the structures, that the resulting orientation should follow
                                     const coord_list &arg_coords         ///< The (pre-superposed) coordinates for which the orientation should be calculated
                                     ) {
	const auto orienting_transformation = get_orienting_transformation( arg_coords );
	post_translate_and_rotate(
		arg_superposition,
		orienting_transformation.first,
		orienting_transformation.second
	);
}

/// \brief Return a copy the specified superposition, modified so as to orient the specified (pre-superposed) coordinates
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
///
/// Note: the order of the coordinates matters
superposition cath::sup::orient_superposition_copy(superposition     arg_superposition, ///< The superposition of the structures, that the resulting orientation should follow
                                                   const coord_list &arg_coords         ///< The (pre-superposed) coordinates for which the orientation should be calculated
                                                   ) {
	orient_superposition( arg_superposition, arg_coords );
	return arg_superposition;
}

/// \brief Modify the specified superposition so as to superpose the specified alignment/PDBs
///        as before but with the whole superposition "correctly" post-oriented
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
void cath::sup::orient_superposition(superposition   &arg_superposition, ///< The superposition of the structures, on top of which the computed orientation should be applied
                                     const alignment &arg_alignment,     ///< The alignment defining the coordinates (including their order) to be oriented
                                     const pdb_list  &arg_pdbs           ///< The PDBs whose coordinates (after superposing) are to be oriented
                                     ) {
	const auto orienting_transformation = get_orienting_transformation(
		arg_superposition,
		arg_alignment,
		arg_pdbs
	);
	post_translate_and_rotate(
		arg_superposition,
		orienting_transformation.first,
		orienting_transformation.second
	);
}

/// \brief Return a copy the specified superposition, modified so as to superpose the specified alignment/PDBs
///        as before but with the whole superposition "correctly" post-oriented
///
/// These orient functions aim to put the primary PCA direction on the x-axis and the secondary on the y-axis;
/// they then choose whether to apply further x, y or z flips by seeing which (if any) best pulls the later coords to
/// the positive x/y quadrant (which makes the result independent of the structures' initial orientations and more stable).
superposition cath::sup::orient_superposition_copy(superposition    arg_superposition, ///< The superposition of the structures, on top of which the computed orientation should be applied
                                                   const alignment &arg_alignment,     ///< The alignment defining the coordinates (including their order) to be oriented
                                                   const pdb_list  &arg_pdbs           ///< The PDBs whose coordinates (after superposing) are to be oriented
                                                   ) {
	orient_superposition( arg_superposition, arg_alignment, arg_pdbs );
	return arg_superposition;
}
