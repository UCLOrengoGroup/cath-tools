/// \file
/// \brief The dssp_file class definitions

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

#include "dssp_file.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/count_if.hpp>

#include "biocore/residue_id.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "file/dssp_wolf/tally_residue_ids.hpp"
#include "file/pdb/dssp_skip_policy.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "file/pdb/protein_info.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::lexical_cast;
using boost::numeric_cast;
using boost::range::count_if;

/// \brief Ctor for dssp_file
dssp_file::dssp_file(residue_vec arg_dssp_residues ///< TODOCUMENT
                     ) : dssp_residues { std::move( arg_dssp_residues ) } {
}

/// \brief TODOCUMENT
size_t dssp_file::get_num_residues() const {
	return dssp_residues.size();
}

/// \brief TODOCUMENT
const residue & dssp_file::get_residue_of_index(const size_t &arg_index ///< TODOCUMENT
                                                ) const {
	return dssp_residues[arg_index];
}

/// \brief TODOCUMENT
dssp_file::const_iterator dssp_file::begin() const {
	return cath::common::cbegin( dssp_residues );
}

/// \brief TODOCUMENT
dssp_file::const_iterator dssp_file::end() const {
	return cath::common::cend( dssp_residues );
}

/// \brief Combine a dssp_file and pdb representing the same structure in a sensible protein object
///
/// \relates dssp_file
///
/// \TODO Consider taking an ostream_ref_opt argument rather than assuming cerr
///       (fix all errors, *then* provide default of boost::none)
protein cath::file::protein_from_dssp_and_pdb(const dssp_file        &arg_dssp_file,        ///< The dssp_file object for a given structure
                                              const pdb              &arg_pdb_file,         ///< The dssp_file object for a given structure
                                              const dssp_skip_policy &arg_dssp_skip_policy, ///< Whether to exclude residues that are in the PDB but not the DSSP
                                              const string           &arg_name,             ///< The name to set as the title of the protein
                                              const ostream_ref_opt  &arg_ostream           ///< An optional reference to an ostream to which any logging should be sent
                                              ) {
	// Build a rough protein object from the pdb object
	const auto  built_protein     = build_protein_of_pdb(
		arg_pdb_file,
		arg_ostream,
		( arg_dssp_skip_policy == dssp_skip_policy::SKIP__BREAK_ANGLES )
			? dssp_skip_policy::DONT_SKIP__BREAK_ANGLES
			: arg_dssp_skip_policy
	);
	const auto &pdb_protein       = built_protein.first;
	const auto &pdb_prot_info     = built_protein.second;
	const auto  pdb_skip_indices  = get_protein_res_indices_that_dssp_might_skip( arg_pdb_file, arg_ostream );

	// Grab the number of residues in the protein and dssp_file objects
	const auto  num_dssp_residues = arg_dssp_file.get_num_residues();
	const auto  num_pdb_residues  = pdb_protein.get_length();

	// Grab the residues names from the DSSP and PDB and then tally them up
	const auto  pdb_res_names     = get_residue_ids  ( pdb_protein );
	const auto  dssp_res_names    = get_residue_ids  ( arg_dssp_file, false );
	const auto  alignment         = tally_residue_ids(
		pdb_res_names,
		dssp_res_names,
		false,
		true,
		pdb_skip_indices
	);

	// Prepare a list of new residue to populate
	residue_vec new_residues;
	new_residues.reserve( ( arg_dssp_skip_policy == dssp_skip_policy::SKIP__BREAK_ANGLES ) ? num_dssp_residues : num_pdb_residues );

	// Loop over the residues
	size_t alignment_ctr = 0;
	for (const size_t &pdb_residue_ctr : indices( num_pdb_residues ) ) {
		const residue        &the_pdb_residue = pdb_protein.get_residue_ref_of_index( pdb_residue_ctr );
		const residue_makeup &res_makeup      = pdb_prot_info.residue_makeups[ pdb_residue_ctr ];

		// If this PDB residue is in the alignment then it can be combined with the equivalent DSSP residue
		const bool is_in_alignment     = ( (alignment_ctr < alignment.size() ) && ( alignment[alignment_ctr].first == pdb_residue_ctr ) );
		if ( is_in_alignment ) {
			// Combine the two residues and add them to the back
			const residue &the_dssp_residue = arg_dssp_file.get_residue_of_index( alignment[alignment_ctr].second );
			new_residues.push_back(
				combine_residues_from_dssp_and_pdb(
					the_dssp_residue,
					the_pdb_residue,
					angle_skipping_of_dssp_skip_policy( arg_dssp_skip_policy ),
					res_makeup
				)
			);

			// Increment the alignment counter
			++alignment_ctr;
		}
		else if ( res_skipping_of_dssp_skip_policy( arg_dssp_skip_policy ) == dssp_skip_res_skipping::DONT_SKIP ) {
			new_residues.push_back( the_pdb_residue );
		}
	}

	// Construct a new protein from the new list of residues
	return { arg_name, new_residues };
}

/// \brief Extract the names of the residues of the specified dssp_file
///
/// \relates dssp_file
///
/// \returns A vector of strings, each containing the PDB residue name or an empty string for a null residue record
residue_id_vec cath::file::get_residue_ids(const dssp_file &arg_dssp_file,            ///< The dssp_file object for a given structure
                                           const bool      &arg_exclude_null_residues ///< Whether to exclude null residues (rather than represent them with empty strings)
                                           ) {
	return copy_build<residue_id_vec>(
		arg_dssp_file
			// Include if this isn't a null residue or if null residues aren't to be excluded
			| filtered(
				[&] (const residue &x) {
					return ! arg_exclude_null_residues || ! is_null_residue( x );
				}
			)
			// Grab the residue name of the residue (or residue_name() for null residues)
			| transformed(
				[] (const residue &x) {
					return is_null_residue( x ) ? residue_id{} : x.get_pdb_residue_id();
				}
			)
	);
}

/// \brief TODOCUMENT
///
/// \relates dssp_file
size_t cath::file::get_num_non_null_residues(const dssp_file &arg_dssp_file ///< TODOCUMENT
                                             ) {
	return numeric_cast<size_t>( count_if(
		arg_dssp_file,
		[] (const residue &x) { return ! is_null_residue( x ); }
	) );
}
