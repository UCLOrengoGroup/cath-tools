/// \file
/// \brief The alignment_coord_extractor class definitions

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

#include "alignment_coord_extractor.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>

#include "alignment/align_type_aliases.hpp"
#include "alignment/common_atom_selection_policy/common_atom_selection_policy.hpp"
#include "alignment/common_residue_selection_policy/common_residue_selection_policy.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/protein_info.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using std::cref;

alignment_coord_extractor::residue_cref_residue_cref_pair_vec alignment_coord_extractor::get_common_residues(const alignment                       &prm_alignment,        ///< The alignment to determine which residues should be as close as possible to which
                                                                                                             const protein                         &prm_protein_a,        ///< Coordinates for first structure
                                                                                                             const protein                         &prm_protein_b,        ///< Coordinates for second structure
                                                                                                             const common_residue_selection_policy &prm_res_seln_policy,  ///< TODOCUMENT
                                                                                                             const alignment::size_type            &prm_entry_index_a,    ///< TODOCUMENT
                                                                                                             const alignment::size_type            &prm_entry_index_b     ///< TODOCUMENT
                                                                                                             ) {
	using aln_size_type = alignment::size_type;
//	const aln_size_type alignment_length = prm_alignment.length();

	// Use the common_residue_selection_policy to determine which indices to grab
	const vector<aln_size_type> selection_aln_indices = prm_res_seln_policy.select_common_residues(prm_alignment, prm_entry_index_a, prm_entry_index_b);

	// TODOCUMENT
	for (const aln_size_type &alignment_index : selection_aln_indices) {
		const bool aln_has_posn_for_first_entry(  has_position_of_entry_of_index( prm_alignment, prm_entry_index_a, alignment_index) );
		const bool aln_has_posn_for_second_entry( has_position_of_entry_of_index( prm_alignment, prm_entry_index_b, alignment_index) );
		if ( ! aln_has_posn_for_first_entry || ! aln_has_posn_for_second_entry ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to retrieve both entries for alignment index selected by common_residue_selection_policy"));
		}
	}

	// Fill array with coordinates
	residue_cref_residue_cref_pair_vec residues;
	for (const aln_size_type &alignment_index : selection_aln_indices) {
		const aln_posn_type a_position = get_position_of_entry_of_index( prm_alignment, prm_entry_index_a, alignment_index );
		const aln_posn_type b_position = get_position_of_entry_of_index( prm_alignment, prm_entry_index_b, alignment_index );
		const residue      &residue_a  = prm_protein_a.get_residue_ref_of_index( a_position );
		const residue      &residue_b  = prm_protein_b.get_residue_ref_of_index( b_position );

//		if ( ! residue_a.get_pdb_name_number() &&  ! residue_b.get_pdb_name_number() ) {
//			BOOST_THROW_EXCEPTION(runtime_error_exception(
//				"When getting common coordinates from alignment, unable to get pdb_name_number() from either residue selected by policy at alignment position "
//				+ lexical_cast<string>( alignment_index )
//				+ " (residue a position: "
//				+ lexical_cast<string>( a_position      )
//				+ "; residue b position: "
//				+ lexical_cast<string>( a_position      )
//				+ ")"
//			));
//		}
		residues.push_back( make_pair( std::cref( residue_a ), std::cref( residue_b ) ) );
	}

	return residues;
}

/// \brief A static method to get the common coordinates, as defined by an alignment, from two protein objects.
///
/// See alignment_coord_extractor for more information.
pair<coord_list_vec, coord_list_vec> alignment_coord_extractor::get_common_coords_by_residue(const alignment                       &prm_alignment,        ///< The alignment to determine which residues should be as close as possible to which
                                                                                             const protein                         &prm_protein_a,        ///< Coordinates for first structure
                                                                                             const protein                         &prm_protein_b,        ///< Coordinates for second structure
                                                                                             const common_residue_selection_policy &prm_res_seln_policy,  ///< TODOCUMENT
                                                                                             const common_atom_selection_policy    &prm_atom_seln_policy, ///< TODOCUMENT
                                                                                             const alignment::size_type            &prm_entry_index_a,    ///< TODOCUMENT
                                                                                             const alignment::size_type            &prm_entry_index_b     ///< TODOCUMENT
                                                                                             ) {
	const residue_cref_residue_cref_pair_vec residues = get_common_residues(
		prm_alignment,
		prm_protein_a,
		prm_protein_b,
		prm_res_seln_policy,
		prm_entry_index_a,
		prm_entry_index_b
	);

	// Declare 2 vectors of coordinates
	pair<coord_list_vec, coord_list_vec> coord_lists;
	coord_lists.first.reserve(  prm_alignment.length() );
	coord_lists.second.reserve( prm_alignment.length() );

	// Fill array with coordinates
	for (const residue_cref_residue_cref_pair &residue_pair : residues) {
		// Use the common_atom_selection_policy to perform the common atom selections
		const pair<coord_list, coord_list> new_residue_coord_lists = select_common_atoms(
			prm_atom_seln_policy,
			residue_pair.first.get(),
			residue_pair.second.get()
		);
		coord_lists.first.push_back(  new_residue_coord_lists.first  );
		coord_lists.second.push_back( new_residue_coord_lists.second );
	}

	return coord_lists;
}

/// \brief A static method to get the common coordinates, as defined by an alignment, from two protein objects.
///
/// See alignment_coord_extractor for more information.
pair<coord_list, coord_list> alignment_coord_extractor::get_common_coords(const alignment                       &prm_alignment,        ///< The alignment to determine which residues should be as close as possible to which
                                                                          const protein                         &prm_protein_a,        ///< Coordinates for first structure
                                                                          const protein                         &prm_protein_b,        ///< Coordinates for second structure
                                                                          const common_residue_selection_policy &prm_res_seln_policy,  ///< TODOCUMENT
                                                                          const common_atom_selection_policy    &prm_atom_seln_policy, ///< TODOCUMENT
                                                                          const alignment::size_type            &prm_entry_index_a,    ///< TODOCUMENT
                                                                          const alignment::size_type            &prm_entry_index_b     ///< TODOCUMENT
                                                                          ) {
	const residue_cref_residue_cref_pair_vec residues = get_common_residues(
		prm_alignment,
		prm_protein_a,
		prm_protein_b,
		prm_res_seln_policy,
		prm_entry_index_a,
		prm_entry_index_b
	);

	// Declare 2 vectors of coordinates
	pair<coord_list, coord_list> coord_lists;
	coord_lists.first.reserve(  prm_alignment.length() );
	coord_lists.second.reserve( prm_alignment.length() );

	// Fill array with coordinates
	for (const residue_cref_residue_cref_pair &residue_pair : residues) {
		// Use the common_atom_selection_policy to perform the common atom selections
		prm_atom_seln_policy.append_common_atoms_to_coord_lists( coord_lists, residue_pair.first.get(), residue_pair.second.get() );
	}

	return coord_lists;
}

/// \brief A static method to get the common coordinates, as defined by an alignment, from two pdb objects.
///
/// See alignment_coord_extractor for more information
///
/// \TODO Consider taking an ostream_ref_opt argument rather than assuming cerr
///       (fix all errors, *then* provide default of boost::none)
pair<coord_list, coord_list> alignment_coord_extractor::get_common_coords(const alignment                       &prm_alignment,        ///< The alignment to determine which residues should be as close as possible to which
                                                                          const pdb                             &prm_pdb_a,            ///< Coordinates for first structure
                                                                          const pdb                             &prm_pdb_b,            ///< Coordinates for second structure
                                                                          const common_residue_selection_policy &prm_res_seln_policy,  ///< TODOCUMENT
                                                                          const common_atom_selection_policy    &prm_atom_seln_policy, ///< TODOCUMENT
                                                                          const alignment::size_type            &prm_entry_index_a,    ///< TODOCUMENT
                                                                          const alignment::size_type            &prm_entry_index_b     ///< TODOCUMENT
                                                                          ) {
	return get_common_coords(
		prm_alignment,
		build_protein_of_pdb( prm_pdb_a, ref( cerr ) ).first,
		build_protein_of_pdb( prm_pdb_b, ref( cerr ) ).first,
		prm_res_seln_policy,
		prm_atom_seln_policy,
		prm_entry_index_a,
		prm_entry_index_b
	);
}
