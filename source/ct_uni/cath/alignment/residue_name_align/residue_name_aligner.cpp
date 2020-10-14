/// \file
/// \brief The residue_name_aligner class definitions

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

#include "residue_name_aligner.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/alignment/pair_alignment.hpp"
#include "cath/alignment/residue_name_align/detail/residue_name_align_map.hpp"
#include "cath/alignment/residue_score/residue_scorer.hpp"
#include "cath/biocore/residue_name.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <string>

using namespace cath::align;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::none;

/// \brief Construct an alignment between multiple lists of residues by pulling together residues of the same name.
///
/// Does this code necessarily need to construct an alignment
/// (as opposed to just identifying the coordinates to superpose)?
///
/// Advantages:
///  - It helps to illustrate how things are being superposed
///  - It makes it easier to check for (and disallow) conflicts in the ordering of residues
///  - It makes it possible to SSAP-score the alignment and then perform a weighted superposition
///
/// Disadvantages:
///  - It makes things a bit more complicated
///  - Any decisions that have to be taken about how to align the non-matching residues are arbitrary
///    (but may misinterpreted as meaningful).
///
/// Where there are choices about which entry's position to add first:
///  - add the longest entry's position first (helps to give consistent answers for testing)
///  - add the first entry's position first
alignment residue_name_aligner::residue_name_align(const residue_name_vec_vec &prm_residue_lists ///< TODOCUMENT
                                                   ) {
	// Check that there is at least one list
	const str_vec_vec::size_type num_lists = prm_residue_lists.size();
	if (num_lists < 1) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot residue_name_align() zero residue lists"));
	}
	// Check that at least one of the lists is not empty
	// \todo Change this to use C++11 any_of() and use a lambda to check for not empty
	bool found_non_empty = false;
	for (const residue_name_vec &prm_residue_list : prm_residue_lists) {
		if ( ! prm_residue_list.empty() ) {
			found_non_empty = true;
		}
	}
	if ( ! found_non_empty ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot residue_name_align() residue lists that are all empty"));
	}

	// Build a vector of residue_name_align_map objects, one for each prm_residue_list
	vector<residue_name_align_map> maps;
	maps.reserve(num_lists);
	for (const residue_name_vec &residue_list : prm_residue_lists) {
		maps.push_back( make_residue_name_align_map( residue_list ) );
	}

	// Data structures to create the alignment:
	size_vec next_index_to_add_for_lists(num_lists, 0);
	aln_posn_opt_vec_vec raw_alignment_data( num_lists );

	// Do the actual work of building an alignment
	bool more_to_do = true;
	while (more_to_do) {
		more_to_do = false;

		// This is made a bit more complicated to ensure that the longer entry's position
		// is inserted first where there is a choice
		bool found_non_skipping = false;
		size_t entry_of_non_skipping = 0;
		// Search to find an entry for which a new entry can be inserted into the alignment
		for (const size_t &entry_ctr : indices( num_lists ) ) {

			// If this entry is complete, then continue to the next entry
			// otherwise, grab the next residue string for this entry
			// and record that there's more_to_do
			const size_t &next_index_to_add = next_index_to_add_for_lists[entry_ctr];
			if ( next_index_to_add >= prm_residue_lists[ entry_ctr ].size() ) {
				continue;
			}
			const residue_name &the_res_name = prm_residue_lists[entry_ctr][next_index_to_add];
			more_to_do = true;

			// If a better or equal non-skipping has already been found (ie for an entry that's at least as long as this one)
			// then no point considering this one so continue to next pass of loop
			if (found_non_skipping && prm_residue_lists[entry_of_non_skipping].size() >= prm_residue_lists[entry_ctr].size()) {
				continue;
			}

			// Find which entries have equivalent residues and the indices of those equivalents
			aln_posn_opt_vec equivalent_indices;
			equivalent_indices.reserve( num_lists );
			for (const residue_name_align_map &map : maps) {
				const aln_posn_opt value = contains_residue_name( map, the_res_name ) ? aln_posn_opt( get_index_of_residue_name( map, the_res_name ) )
				                                                                      : aln_posn_opt( none );
				equivalent_indices.push_back( value );
			}

			// Check whether inserting this row of equivalents would involve skipping anything
			bool found_skip_here = false;
			for (const size_t &entry_check_ctr : indices( num_lists ) ) {
				const aln_posn_opt &position = equivalent_indices[entry_check_ctr];
				if ( position && *position != next_index_to_add_for_lists[ entry_check_ctr ] ) {
					if ( *position < next_index_to_add_for_lists[entry_check_ctr]) {
						BOOST_THROW_EXCEPTION(
							invalid_argument_exception(
								"Whilst aligning residue names, residue " + lexical_cast<string>( the_res_name    )
								+ " is out of order (in entry "           + lexical_cast<string>( entry_ctr       )
								+ " and entry "                           + lexical_cast<string>( entry_check_ctr )
								+ " around indices "                      + lexical_cast<string>( *position       )
								+ " and "                                 + lexical_cast<string>( next_index_to_add_for_lists[ entry_check_ctr ] )
								+ ")"
							)
						);
					}
					found_skip_here = true;
					break;
				}
			}

			// If we've found a skip here then it's better than before If we've found a better skip than before then record it
			if ( ! found_skip_here ) {
				entry_of_non_skipping = entry_ctr;
				found_non_skipping = true;
			}
			else {
			}
		}

		// If this entry is ready to be inserted, then proceed and break out of this loop
		if ( found_non_skipping ) {
			const residue_name &the_res_name = prm_residue_lists[entry_of_non_skipping][next_index_to_add_for_lists[entry_of_non_skipping]];

			// Find which entries have equivalent residues
			bool_deq equivalent_presences;
			for (const residue_name_align_map &map : maps) {
				equivalent_presences.push_back( contains_residue_name( map, the_res_name ) );
			}

			// Insert the new positions and increment the relevant indices in next_index_to_add_for_lists
			for (const size_t &entry_check_ctr : indices( num_lists ) ) {
				const bool &should_insert_entry = equivalent_presences[ entry_check_ctr ];
				const aln_posn_opt value = should_insert_entry ? next_index_to_add_for_lists[ entry_check_ctr ]
				                                               : aln_posn_opt( none );
				raw_alignment_data[ entry_check_ctr ].push_back( value );
				if ( should_insert_entry ) {
					++( next_index_to_add_for_lists[ entry_check_ctr ] );
				}
			}
		}

		if ( more_to_do && ! found_non_skipping ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to align residue names, this probably means the residue ordering is not consistent"));
		}
	}

	// Create and return an alignment from the data
	return alignment( raw_alignment_data );
}

/// \brief TODOCUMENT
alignment cath::align::residue_name_align_and_residue_score(const residue_name_vec_vec &prm_residue_lists,  ///< TODOCUMENT
                                                            const residue_scorer       &prm_residue_scorer, ///< TODOCUMENT
                                                            const protein_list         &prm_proteins        ///< TODOCUMENT
                                                            ) {
	alignment new_alignment = residue_name_aligner::residue_name_align( prm_residue_lists );
	score_alignment( prm_residue_scorer, new_alignment, prm_proteins );
	return new_alignment;
}

/// \brief TODOCUMENT
alignment cath::align::residue_name_align_and_residue_score_if_multi(const residue_name_vec_vec &prm_residue_lists,  ///< TODOCUMENT
                                                                     const residue_scorer       &prm_residue_scorer, ///< TODOCUMENT
                                                                     const protein_list         &prm_proteins        ///< TODOCUMENT
                                                                     ) {
	const size_t num_proteins = prm_proteins.size();
	return (num_proteins > 1) ? residue_name_align_and_residue_score( prm_residue_lists, prm_residue_scorer, prm_proteins )
	                          : residue_name_aligner::residue_name_align( prm_residue_lists );
}
