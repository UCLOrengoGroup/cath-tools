/// \brief The protein_io class definitions


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
#include "protein_io.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/lexical_cast.hpp>

#include "common/file/open_fstream.hpp"
#include "common/lexical_cast_line.hpp"
#include "common/type_aliases.hpp"
#include "exception/runtime_error_exception.hpp"
#include "file/dssp_wolf/dssp_file.hpp"
#include "file/dssp_wolf/dssp_file_io.hpp"
#include "file/dssp_wolf/wolf_file.hpp"
#include "file/dssp_wolf/wolf_file_io.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "file/sec/sec_file.hpp"
#include "file/sec/sec_file_io.hpp"
#include "file/sec/sec_file_record.hpp"
#include "ssap/clique.hpp"
#include "ssap/ssap.hpp"
#include "structure/accessibility_calc/dssp_accessibility.hpp"
#include "structure/entry_querier/entry_querier.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "structure/sec_struc_calc/dssp/dssp_ss_calc.hpp"
#include "structure/sec_struc_calc/sec/sec_calc.hpp"

#include <iostream>
#include <string>
#include <vector>

using namespace boost::filesystem;
using namespace cath;
using namespace cath::file;
using namespace cath::sec;
using namespace std;

using boost::algorithm::any_of;
using boost::lexical_cast;
using boost::numeric_cast;

/// \brief Read a wolf and a sec file and build them into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_wolf_and_sec_files(const path            &arg_wolf_filename, ///< The WOLF file to read
                                                   const path            &arg_sec_filename,  ///< A sec file
                                                   const string          &arg_name,          ///< The name to set as the title of the protein
                                                   const ostream_ref_opt &arg_stderr         ///< An optional reference to an ostream to which any logging should be performed
                                                   ) {
	return protein_from_wolf_and_sec(
		read_wolf( arg_wolf_filename ),
		read_sec(  arg_sec_filename  ),
		arg_name,
		arg_stderr
	);
}

/// \brief Read a DSSP, a PDB and a sec file and build them into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_dssp_pdb_and_sec_files(const path             &arg_dssp_filename,    ///< A DSSP file
                                                       const path             &arg_pdb_filename,     ///< A PDB file
                                                       const path             &arg_sec_filename,     ///< A sec file
                                                       const dssp_skip_policy &arg_dssp_skip_policy, ///< Whether to limit the protein to those residues that were present in the DSSP
                                                       const string           &arg_name,             ///< The name to set as the title of the protein
                                                       const ostream_ref_opt  &arg_stderr            ///< An optional reference to an ostream to which any logging should be performed
                                                       ) {
	return add_name_and_paint_sec_file_onto_protein_copy(
		read_protein_from_dssp_and_pdb(
			arg_dssp_filename,
			arg_pdb_filename,
			arg_dssp_skip_policy,
			""s,
			arg_stderr
		),
		read_sec(
			arg_sec_filename
		),
		arg_name,
		arg_stderr
	);
}


/// \brief Read a DSSP, a PDB, calculate the corresponding sec_file and build it all into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_dssp_and_pdb_and_calc_sec(const path             &arg_dssp_filename,    ///< A DSSP file
                                                          const path             &arg_pdb_filename,     ///< A PDB file
                                                          const dssp_skip_policy &arg_dssp_skip_policy, ///< Whether to limit the protein to those residues that were present in the DSSP
                                                          const string           &arg_name,             ///< The name to set as the title of the protein
                                                          const ostream_ref_opt  &arg_stderr            ///< An optional reference to an ostream to which any logging should be performed
                                                          ) {
	protein the_protein = read_protein_from_dssp_and_pdb(
		arg_dssp_filename,
		arg_pdb_filename,
		arg_dssp_skip_policy,
		""s,
		arg_stderr
	);
	add_name_and_paint_sec_file_onto_protein(
		the_protein,
		get_sec_file( the_protein ),
		arg_name,
		arg_stderr
	);
	return the_protein;
}

/// \brief Read a DSSP and a PDB file and build them into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_dssp_and_pdb(const path             &arg_dssp,             ///< A DSSP file
                                             const path             &arg_pdb,              ///< A PDB file
                                             const dssp_skip_policy &arg_dssp_skip_policy, ///< Whether to limit the protein to those residues that were present in the DSSP
                                             const string           &arg_name,             ///< The name to set as the title of the protein
                                             const ostream_ref_opt  &arg_ostream           ///< An optional reference to an ostream to which any logging should be sent
                                             ) {
	return protein_from_dssp_and_pdb(
		read_dssp_file( arg_dssp ),
		read_pdb_file ( arg_pdb  ),
		arg_dssp_skip_policy,
		arg_name,
		arg_ostream
	);
}

/// \brief Read a DSSP and a PDB file and build them into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_pdb(const path   &arg_pdb,   ///< A PDB file
                                    const string &arg_name,  ///< The name to set as the title of the protein
                                    ostream      &arg_stderr ///< TODOCUMENT
                                    ) {
	return build_protein_of_pdb_and_name(
		read_pdb_file( arg_pdb ),
		arg_name,
		ref( arg_stderr )
	);
}

/// \brief Read a DSSP and a PDB file and build them into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_pdb_and_calc_dssp(const path            &arg_pdb,    ///< A PDB file
                                                  const string          &arg_name,   ///< The name to set as the title of the protein
                                                  const ostream_ref_opt &arg_ostream ///< An optional reference to an ostream to which any logging should be sent
                                                  ) {
	const pdb the_pdb = backbone_complete_subset_of_pdb(
		read_pdb_file( arg_pdb ),
		arg_ostream
	).first;
	protein the_protein = build_protein_of_pdb_and_name(
		the_pdb,
		arg_name,
		arg_ostream
	);

	// build_protein_of_pdb_and_name() already does phi/psi angles
	// so now do the sec_struc types and accessibilities
	set_sec_struc_types( the_protein, calc_sec_strucs_of_backbone_complete_pdb( the_pdb ) );
	set_accessibilities( the_protein, calc_accessibilities_with_scanning      ( the_pdb ) );
	return the_protein;
}

/// \brief Read a DSSP and a PDB file and build them into a protein
///
/// \relatesalso protein
protein cath::read_protein_from_pdb_and_calc_dssp_and_sec(const path            &arg_pdb,    ///< A PDB file
                                                          const string          &arg_name,   ///< The name to set as the title of the protein
                                                          const ostream_ref_opt &arg_ostream ///< An optional reference to an ostream to which any logging should be sent
                                                          ) {
	protein the_protein = read_protein_from_pdb_and_calc_dssp(
		arg_pdb,
		arg_name,
		arg_ostream
	);
	add_name_and_paint_sec_file_onto_protein(
		the_protein,
		get_sec_file( the_protein ),
		arg_name,
		arg_ostream
	);
	return the_protein;
}

/// \brief Construct a protein object from parsed WOLF and sec files
///
/// \relatesalso protein
/// \relatesalso wolf_file
/// \relatesalso sec_file
protein cath::protein_from_wolf_and_sec(const wolf_file       &arg_wolf,  ///< The parsed wolf_file object
                                        const sec_file        &arg_sec,   ///< The parsed sec
                                        const string          &arg_name,  ///< The name to set as the title of the protein
                                        const ostream_ref_opt &arg_stderr ///< An optional reference to an ostream to which any logging should be performed
                                        ) {
	// Create a protein from the wolf file, set its title, paint the sec_file onto it and then return it
	return add_name_and_paint_sec_file_onto_protein_copy(
		protein_from_wolf(arg_wolf),
		arg_sec,
		arg_name,
		arg_stderr
	);
}

/// \brief Construct a protein object from parsed DSSP, PDB and sec files
///
/// \relatesalso protein
/// \relatesalso dssp_file
/// \relatesalso pdb
/// \relatesalso sec_file
protein cath::protein_from_dssp_pdb_and_sec(const dssp_file        &arg_dssp,             ///< The parsed DSSP file
                                            const pdb              &arg_pdb,              ///< The parsed PDB file
                                            const sec_file         &arg_sec,              ///< The parsed sec file
                                            const dssp_skip_policy &arg_dssp_skip_policy, ///< Whether to limit the protein to those residues that were present in the DSSP
                                            const string           &arg_name,             ///< The name to set as the title of the protein
                                            const ostream_ref_opt  &arg_stderr            ///< An optional reference to an ostream to which any logging should be performed
                                            ) {
	// Create a protein from the DSSP and PDB, set its title, paint the sec_file onto it and then return it
	return add_name_and_paint_sec_file_onto_protein_copy(
		protein_from_dssp_and_pdb(
			arg_dssp,
			arg_pdb,
			arg_dssp_skip_policy,
			arg_name,
			arg_stderr
		),
		arg_sec,
		arg_name,
		arg_stderr
	);
}

/// \brief Paint the sec_file's secondary structures onto the protein and set its name
///
/// \relatesalso protein
/// \relatesalso sec_file
void cath::add_name_and_paint_sec_file_onto_protein(protein               &arg_protein, ///< The protein to be modified (taken by non-const reference to allow it to be modified)
                                                    const sec_file        &arg_sec,     ///< The parsed sec file whose secondary structures should be painted on to the protein
                                                    const string          &arg_name,    ///< The name to set as the title of the protein
                                                    const ostream_ref_opt &arg_stderr   ///< An optional reference to an ostream to which any logging should be performed
                                                    ) {
	// Set the protein's title and then paint the sec_file onto it
	arg_protein.set_title(arg_name);
	paint_sec_file_onto_protein( arg_protein, arg_sec, arg_stderr );
}

/// \brief Paint the sec_file's secondary structures onto a copy of the protein and set its name, then return it
///
/// \relatesalso protein
/// \relatesalso sec_file
protein cath::add_name_and_paint_sec_file_onto_protein_copy(protein                arg_protein, ///< The protein whose copy should be modified (taken by value to allow the compiler to elide copies of rvalues)
                                                            const sec_file        &arg_sec,     ///< The parsed sec file whose secondary structures should be painted on to the protein
                                                            const string          &arg_name,    ///< The name to set as the title of the protein
                                                            const ostream_ref_opt &arg_stderr   ///< An optional reference to an ostream to which any logging should be performed
                                                            ) {
	// Using the local copy of the protein, set it's title and paint the sec_file onto it and then return it
	add_name_and_paint_sec_file_onto_protein( arg_protein, arg_sec, arg_name, arg_stderr );
	return arg_protein;
}

/// \brief Paint the sec_file's secondary structures onto the protein
///
/// \relatesalso protein
/// \relatesalso sec_file
void cath::paint_sec_file_onto_protein(protein               &arg_protein,  ///< The protein to be modified (taken by non-const reference to allow it to be modified)
                                       const sec_file        &arg_sec_file, ///< The parsed sec file whose secondary structures should be painted on to the protein
                                       const ostream_ref_opt &arg_stderr    ///< An optional reference to an ostream to which any logging should be performed
                                       ) {
	const sec_struc_vec new_sec_strucs = make_sec_struc_list(arg_sec_file);
	arg_protein.set_sec_strucs ( new_sec_strucs );

	// Label residues with secondary structure number
	label_residues_with_sec_strucs(arg_protein, arg_stderr);
}

/// \brief Paint the sec_file's secondary structures onto a copy of the protein, then return it
///
/// \relatesalso protein
/// \relatesalso sec_file
protein cath::paint_sec_file_onto_protein_copy(protein                arg_protein,  ///< The protein whose copy should be modified (taken by value to allow the compiler to elide copies of rvalues)
                                               const sec_file        &arg_sec_file, ///< The parsed sec file whose secondary structures should be painted on to the protein
                                               const ostream_ref_opt &arg_stderr    ///< An optional reference to an ostream to which any logging should be performed
                                               ) {
	paint_sec_file_onto_protein(arg_protein, arg_sec_file, arg_stderr);
	return arg_protein;
}

/// \brief Calculate a sec_file from the specified protein and paint its secondary structures onto the protein
///
/// \relatesalso protein
/// \relatesalso sec_file
void cath::calc_and_paint_sec_file_onto_protein(protein               &arg_protein, ///< The protein to be modified (taken by non-const reference to allow it to be modified)
                                                const ostream_ref_opt &arg_stderr   ///< An optional reference to an ostream to which any logging should be performed
                                                ) {
	arg_protein.set_sec_strucs( make_sec_struc_list( get_sec_file( arg_protein ) ) );

	// Label residues with secondary structure number
	label_residues_with_sec_strucs( arg_protein, arg_stderr );
}

/// \brief Calculate a sec_file from the specified protein and paint its secondary structures onto a copy of the protein, the return it
///
/// \relatesalso protein
/// \relatesalso sec_file
protein cath::calc_and_paint_sec_file_onto_protein_copy(protein                arg_protein, ///< The protein whose copy should be modified (taken by value to allow the compiler to elide copies of rvalues)
                                                        const ostream_ref_opt &arg_stderr   ///< An optional reference to an ostream to which any logging should be performed
                                                        ) {
	calc_and_paint_sec_file_onto_protein( arg_protein, arg_stderr);
	return arg_protein;
}

/// \brief Remove domin regions to blank out of search chain
///
/// \relatesalso protein
void cath::remove_domin_res(protein               &arg_protein,        ///< The protein object to modify
                            const path            &arg_domin_filename, ///< A domin file
                            const ostream_ref_opt &arg_stderr          ///< An optional reference to an ostream to which any logging should be performed
                            ) {
	// Prepare a vector to contain the starts and ends from the clique data
	size_size_pair_vec clique_starts_and_ends;

	// Read the domin file
	clique new_clique_file( read_clique_file( arg_domin_filename ) );

	// Grab the data of interest from the clique object
	const size_t clique_size = new_clique_file.cliquesize;
	clique_starts_and_ends.reserve(clique_size);
	for (size_t clique_ctr = 0; clique_ctr < clique_size; ++clique_ctr) {
		clique_starts_and_ends.push_back(
			make_pair(
				lexical_cast<size_t>( new_clique_file.equivs[ clique_ctr ].prota_start ),
				lexical_cast<size_t>( new_clique_file.equivs[ clique_ctr ].prota_end   )
			)
		);
	}

	remove_domin_res( arg_protein, clique_starts_and_ends, arg_stderr );
}

/// \brief Remove domin regions to blank out of search chain
///
/// \relatesalso protein
void cath::remove_domin_res(protein                  &arg_protein,                ///< The protein object to modify
                            const size_size_pair_vec &arg_clique_starts_and_ends, ///< The starts and stops that have been parsed out of a domin file
                            const ostream_ref_opt    &arg_stderr                  ///< An optional reference to an ostream to which any logging should be performed
                            ) {
	/////
	// 1. Process the secondary structures that are to be removed
	/////

	const size_t num_sec_strucs = arg_protein.get_num_sec_strucs();

	// A vector to populate with the sec_strucs that are to be kept
	sec_struc_vec sec_strucs_to_keep;
	sec_strucs_to_keep.reserve(num_sec_strucs);

	// A vector of the old numbers of each of the kept sec_strucs indexed by their new numbers.
	vector<size_t> sec_struc_index_conv;
	sec_struc_index_conv.reserve(num_sec_strucs);

	// Look for secondary structures to exclude
	for (size_t sec_struc_ctr = 0; sec_struc_ctr < num_sec_strucs; ++sec_struc_ctr) {
		const sec_struc &my_sec_struc = arg_protein.get_sec_struc_ref_of_index(sec_struc_ctr);

		// Look to see if we want to exclude this sec_struc
		const bool found = any_of(
			arg_clique_starts_and_ends,
			[&] (const size_size_pair &x) {
				const size_t &a_start = x.first;
				const size_t &a_end   = x.second;

				// If the secondary structure's start or end is within the range, then mark it at as one to exclude
				// (Should this code also exclude secondary structures that straddle the range?
				//  ie (my_sec_struc.from <= a_end) || (my_sec_struc.to >= a_start) ? )
				return (
					( my_sec_struc.get_start_residue_num() >= a_start && my_sec_struc.get_start_residue_num() <= a_end )
					||
					( my_sec_struc.get_stop_residue_num()  >= a_start && my_sec_struc.get_stop_residue_num()  <= a_end ) );
			}
		);

		// If it's not designated to be removed, then keep
		if ( ! found ) {
			// Record the old sec_struc number according to the new sec_struc number
			// (the new sec_struc number is implied because sec_struc_index_conv is populated synchronously with sec_strucs_to_keep)
			sec_struc_index_conv.push_back(sec_struc_ctr);

			// Add the new sec_struc to sec_strucs_to_keep, first updating its sec_struc number to the new value
			sec_struc new_sec_struc(my_sec_struc);
			sec_strucs_to_keep.push_back(my_sec_struc);
		}
	}

	// Update pair data
	const size_t num_sec_strucs_to_keep = sec_strucs_to_keep.size();
	for (size_t new_sec_struc_ctr_i = 0; new_sec_struc_ctr_i < num_sec_strucs_to_keep; ++new_sec_struc_ctr_i) {
		const size_t    &old_sec_struc_index_i = sec_struc_index_conv[new_sec_struc_ctr_i];
		const sec_struc &old_sec_struc         = arg_protein.get_sec_struc_ref_of_index( old_sec_struc_index_i );

		sec_struc_planar_angles_vec new_planar_angles;
		new_planar_angles.reserve(num_sec_strucs_to_keep);

		for (const size_t &old_sec_struc_index_j : sec_struc_index_conv) {
			const sec_struc_planar_angles old_planar_angles = old_sec_struc.get_planar_angles_of_index( old_sec_struc_index_j );

			// The previous code appeared to be wiping all planar_angles here!
			// Not sure what was going on there.
			//                                        Tony Lewis, 30th September 2012
			new_planar_angles.push_back(old_planar_angles);
		}

		sec_struc &new_sec_struc = sec_strucs_to_keep[new_sec_struc_ctr_i];
		new_sec_struc.set_planar_angles(new_planar_angles);
	}
	arg_protein.set_sec_strucs(sec_strucs_to_keep);



	/////
	// 2. Process the residues that are to be removed
	/////

	const size_t protein_length = arg_protein.get_length();
	residue_vec residues_to_keep;
	residues_to_keep.reserve(protein_length);

	// Look for residues to exclude
	for (size_t residue_ctr = 0; residue_ctr < protein_length; ++residue_ctr) {
		const residue &my_residue = arg_protein.get_residue_ref_of_index(residue_ctr);
		bool           found      = false;

		// Look to see if we want to exclude this residue
		for (const size_size_pair &clique_start_and_end : arg_clique_starts_and_ends) {
			const size_t &a_start = clique_start_and_end.first;
			const size_t &a_end   = clique_start_and_end.second;

			/// Problem in which residues names are compared numerically
			/// --------------------------------------------------------
			///
			/// \todo Fix this! this is comparing residue names without handling insert-codes and without any consideration
			/// that the numbers in residue names can be completely out of order.
			///
			/// Scanning a chain with lots of insert-codes, eg:
			///
			///     > cathedral_scan -d 2qriB -c /cath/data/v3_5_0/release_data/Class1-4.s35.cathedral.library --clique --dir /cath/data/v3_5_0/grath
			///
			/// ...produces files like this:
			///
			///     > cat 2qriB1nepA00.clique
			///     5
			///      2   21   30  1    5    8
			///      3   36   41  7   62   65
			///      8   78   83  6   52   59
			///     29  233  236  5   35   44
			///     30  241  250  3   16   21
			///
			/// ...but then CATHEDRAL can hardly be expected to get this right because the grath file is wrong:
			///
			///     > head -n 6 /cath/data/v3_5_0/grath/2qriB.gth
			///     2qriB
			///     32
			///     6 11
			///     21 30
			///     36 41
			///     43 46
			///
			/// Perhaps it'd be better if everything used sequential numbers rather than insert codes.
			///
			/// \todo Is the initial `pdb_number( my_residue ) != 0` check just an error?
			if ( ( pdb_number( my_residue ) != 0 ) && pdb_number( my_residue ) >= numeric_cast<int>( a_start ) && pdb_number( my_residue ) <= numeric_cast<int>( a_end ) ) {
				found = true;
				break;
			}
		}

		// If it's not designated to be removed, then keep
		if (!found) {
			residues_to_keep.push_back(my_residue);
		}
	}

//	const size_t num_residues_to_keep = residues_to_keep.size();
//
//	for (size_t i = 1; i <= num_residues_to_keep; ++i) {
//		residue &newp = residues_to_keep[i-1];
//
//		// Disulfide information (not sure it's important)
//		if ( islower(newp.get_amino_acid()) ) {
//			cerr << "IMPORTANT NOTICE : ******** THIS BIT OF CODE IS NOT COMPLETELY USELESS - HURRAH - PLEASE TELL SOMEONE *******" << endl;
//			newp.set_amino_acid('C');
//		}
//	}

	arg_protein.set_residues( residues_to_keep );

	label_residues_with_sec_strucs( arg_protein, arg_stderr );
}
