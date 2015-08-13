/// \file
/// \brief The superposition_io class definitions

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


#include "superposition_io.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/range/irange.hpp>

#include "common/file/open_fstream.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h"
#include "exception/runtime_error_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"
#include "superposition/superposition.h"

#include <fstream>

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;
using namespace std;

using boost::irange;
using boost::property_tree::json_parser::write_json;
using boost::property_tree::ptree;

/// \brief TODOCUMENT
///
/// \relates superposition
//double cath::get_rmsd_from_pairwise_superposition(const superposition &arg_superposition
//                                                  ) {
//	check_superposition_is_pairwise(arg_superposition);
//	return arg_superposition.get_rmsd_between_index_and_next(superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION);
//}

/// \brief Write a the superposition's essential data out to an XML superposition file
///
/// \relates superposition
///
/// \todo Improve the use of C++. Read:
///        - http://www.ibm.com/developerworks/xml/library/x-ctlbx/index.html
///        - http://www.grinninglizard.com/tinyxml2docs/index.html
///       (tinyxml is available as installable package on orengobuild64)
void cath::sup::write_xml_sup(ostream              &arg_ostream,       ///< TODOCUMENT
                              const superposition  &arg_superposition, ///< TODOCUMENT
                              const str_vec        &arg_ids            ///< TODOCUMENT
                              ) {
	if (arg_superposition.get_num_entries() != superposition::NUM_ENTRIES_IN_PAIRWISE_SUPERPOSITION) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("Currently only able to output XML superposition for two entries"));
	}
	const size_t num_entries = arg_superposition.get_num_entries();
	if (arg_ids.size() != num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of IDs must match the number of entries in the superposition"));
	}

	arg_ostream << setprecision(6) << fixed;
	arg_ostream << "<?xml version=\"1.0\"?>" << endl;
	arg_ostream << "<root>" << endl;

	for (size_t entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		const coord centre_of_gravity = -arg_superposition.get_translation_of_index(entry_ctr);
		arg_ostream << "  <structure" << entry_ctr+1 << " id=\"" << arg_ids[entry_ctr];
		arg_ostream << "\">" << endl;
		arg_ostream << "    <centre x=\"" << centre_of_gravity.get_x();
		arg_ostream <<          "\" y=\"" << centre_of_gravity.get_y();
		arg_ostream <<          "\" z=\"" << centre_of_gravity.get_z();
		arg_ostream <<          "\" />" << endl;
		arg_ostream << "  </structure" << entry_ctr+1 << ">" << endl;
	}

	for (size_t rotation_ctr = 1; rotation_ctr < num_entries; ++rotation_ctr) {
		const rotation sup_rotn = arg_superposition.get_rotation_of_index(rotation_ctr);
		for (size_t dim_ctr_1 = 0; dim_ctr_1 < coord::NUM_DIMS; ++dim_ctr_1) {
			arg_ostream << "  <rotationmatrix name=\"row" << dim_ctr_1+1 << "\"";
			for (size_t dim_ctr_2 = 0; dim_ctr_2 < coord::NUM_DIMS; ++dim_ctr_2) {
				// NOTE: transposing the matrix to be consistent with legacy behaviour
				arg_ostream << " col" << dim_ctr_2+1 << "=\"" << sup_rotn.get_value(dim_ctr_2, dim_ctr_1) << "\"";
			}
			arg_ostream << " />" << endl;
		}
	}
	arg_ostream << "</root>" << endl;
}


/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::write_xml_sup_filename(const superposition  &arg_superposition, ///< TODOCUMENT
                                       const path           &arg_filename,      ///< TODOCUMENT
                                       const str_vec        &arg_ids            ///< TODOCUMENT
                                       ) {
	// A try block with the aim of catching any I/O errors
	try {
		ofstream xml_ostream;
		open_ofstream(xml_ostream, arg_filename);
		write_xml_sup(xml_ostream, arg_superposition, arg_ids);
		xml_ostream.close();
	}
	catch (const boost::exception &ex) {
		throw;
	}
	catch (const std::exception & ex) {
		const string error_message("Cannot write XML superposition file \"" + arg_filename.string() + "\" [" + ex.what() + "] ");
		perror(error_message.c_str());
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_message));
	};
}

///// \brief TODOCUMENT
/////
///// \relates superposition
//ostream & cath::sup::write_superposed_pdb_to_ostream(ostream             &arg_os,            ///< TODOCUMENT
//                                                     const superposition &arg_superposition, ///< TODOCUMENT
//                                                     const pdb           &arg_pdb            ///< TODOCUMENT
//                                                     ) {
//	write_superposed_pdb_to_ostream(
//		arg_os,
//		arg_superposition,
//		arg_pdb,
//		0,
//		false
//	);
//	return arg_os;
//}

/// \brief TODOCUMENT
///
/// \relates superposition
ostream & cath::sup::write_superposed_pdb_to_ostream(ostream             &arg_os,            ///< TODOCUMENT
	                                                 const superposition &arg_superposition, ///< TODOCUMENT
	                                                 pdb                  arg_pdb,           ///< TODOCUMENT
	                                                 const size_t        &arg_chain_index,   ///< TODOCUMENT
	                                                 const bool          &arg_relabel_chain  ///< TODOCUMENT
	                                                 ) {
	// Translate PDB
	arg_pdb += arg_superposition.get_translation_of_index( arg_chain_index );

	// Rotate PDB
	arg_pdb.rotate( arg_superposition.get_rotation_of_index( arg_chain_index ) );

	if ( arg_relabel_chain ) {
		// Label structure with sensible chain label
		//
		// \todo Change this so that it does not relabel chains unless that is required
		//       (when writing a single old-style superposition PDB)
		//
		// *****
		// * THIS IS USED BY THE GENOME3D SUPERPOSITIONS SO ADD A TEST TO
		// * PRESERVE THAT BEHAVIOUR IF ALTERING THIS CODE
		// *****
		if ( arg_chain_index >= superposition::SUPERPOSITION_CHAIN_LABELS.size() ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Ran out of chain labels when relabelling chains"));
		}
		arg_pdb.set_chain_label( superposition::SUPERPOSITION_CHAIN_LABELS[ arg_chain_index ] );
	}

	write_pdb_file(arg_os, arg_pdb);
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates superposition
ostream & cath::sup::write_superposed_pdbs_to_ostream(ostream             &arg_os,            ///< TODOCUMENT
                                                      const superposition &arg_superposition, ///< TODOCUMENT
                                                      const pdb_list       arg_pdbs,          ///< TODOCUMENT
                                                      const bool          &arg_write_script,  ///< TODOCUMENT
                                                      const bool          &arg_relabel_chain  ///< TODOCUMENT
                                                      ) {
	// Sanity check inputs
	const size_t num_entries = arg_superposition.get_num_entries();
	if (num_entries != arg_pdbs.size()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot superpose PDBs because the number of PDBs does not match the number of entries in the superposition"));
	}
	if (num_entries > superposition::SUPERPOSITION_CHAIN_LABELS.size()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot superpose PDBs because there aren't enough chain labels"));
	}

	if (arg_write_script) {
		arg_os << "#!rasmol -script"                        << endl;
		arg_os << "zap"                                     << endl;
		arg_os << "load inline"                             << endl;
		arg_os << "wireframe off"                           << endl;
		arg_os << "select all"                              << endl;
		arg_os << "cartoon"                                 << endl;
		arg_os << "select *"                                << endl;
		arg_os << "color chain"                             << endl;
		arg_os << "exit"                                    << endl;
	}

	// Translate each PDB to midpoint, based on CoG of equivalent positions
	// (This is important for structures with low overlap)
	for (size_t pdb_ctr = 0; pdb_ctr < num_entries; ++pdb_ctr) {
		write_superposed_pdb_to_ostream( arg_os, arg_superposition, arg_pdbs[pdb_ctr], pdb_ctr, arg_relabel_chain );
	}

	return arg_os;
}

///// TODOCUMENT
/////
///// \relates superposition
/////
///// \todo Refactor out common file-opening/exception-handling functionality from these file readers/writers
//void cath::sup::write_superposed_pdb_to_file(const superposition &arg_superposition, ///< TODOCUMENT
//                                             const path          &arg_filename,      ///< TODOCUMENT
//                                             const pdb           &arg_pdb            ///< TODOCUMENT
//                                             ) {
//	ofstream superposition_ostream;
//	open_ofstream(superposition_ostream, arg_filename);
//
//	// Try here to catch any I/O exceptions
//	try {
//		write_superposed_pdb_to_ostream(
//			superposition_ostream,
//			arg_superposition,
//			arg_pdb,
//			0,
//			false
//		);
//
//		// Close the file
//		superposition_ostream.close();
//	}
//	// Catch and immediately rethrow any boost::exceptions
//	// (so that it won't get caught in the next block if it's a std::exception)
//	catch (const boost::exception &ex) {
//		throw;
//	}
//	// Catch any I/O exceptions
//	catch (const std::exception &ex) {
//		const string error_message("Cannot write superposition to file \"" + arg_filename.string() + "\" [" + ex.what() + "] ");
//		perror(error_message.c_str());
//		BOOST_THROW_EXCEPTION(runtime_error_exception(error_message));
//	};
//}

/// TODOCUMENT
///
/// \relates superposition
///
/// \todo Refactor out common file-opening/exception-handling functionality from these file readers/writers
void cath::sup::write_superposed_pdb_to_file(const superposition &arg_superposition, ///< TODOCUMENT
                                             const path          &arg_filename,      ///< TODOCUMENT
                                             const pdb           &arg_pdb,           ///< TODOCUMENT
                                             const size_t        &arg_chain_index,   ///< TODOCUMENT
                                             const bool          &arg_relabel_chain ///< TODOCUMENT
                                             ) {
	ofstream superposition_ostream;
	open_ofstream(superposition_ostream, arg_filename);

	// Try here to catch any I/O exceptions
	try {
		write_superposed_pdb_to_ostream(
			superposition_ostream,
			arg_superposition,
			arg_pdb,
			arg_chain_index,
			arg_relabel_chain
		);

		// Close the file
		superposition_ostream.close();
	}
	// Catch and immediately rethrow any boost::exceptions
	// (so that it won't get caught in the next block if it's a std::exception)
	catch (const boost::exception &ex) {
		throw;
	}
	// Catch any I/O exceptions
	catch (const std::exception &ex) {
		const string error_message("Cannot write superposition to file \"" + arg_filename.string() + "\" [" + ex.what() + "] ");
		perror(error_message.c_str());
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_message));
	};
}

/// TODOCUMENT
///
/// \relates superposition
///
/// \todo Refactor out common file-opening/exception-handling functionality from these file readers/writers
void cath::sup::write_superposed_pdb_to_file(const superposition &arg_superposition, ///< TODOCUMENT
                                             const path          &arg_filename,      ///< TODOCUMENT
                                             const pdb_list      &arg_pdbs,          ///< TODOCUMENT
                                             const bool          &arg_write_script,  ///< TODOCUMENT
                                             const bool          &arg_relabel_chain  ///< TODOCUMENT
                                             ) {
	ofstream superposition_ostream;
	open_ofstream(superposition_ostream, arg_filename);

	// Try here to catch any I/O exceptions
	try {
		write_superposed_pdbs_to_ostream(
			superposition_ostream,
			arg_superposition,
			arg_pdbs,
			arg_write_script,
			arg_relabel_chain
		);

		// Close the file
		superposition_ostream.close();
	}
	// Catch and immediately rethrow any boost::exceptions
	// (so that it won't get caught in the next block if it's a std::exception)
	catch (const boost::exception &ex) {
		throw;
	}
	// Catch any I/O exceptions
	catch (const std::exception &ex) {
		const string error_message("Cannot write superposition to file \"" + arg_filename.string() + "\" [" + ex.what() + "] ");
		perror(error_message.c_str());
		BOOST_THROW_EXCEPTION(runtime_error_exception(error_message));
	};
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::write_superposed_pdb_from_files(const superposition &arg_superposition, ///< TODOCUMENT
                                                const path          &arg_filename,      ///< TODOCUMENT
                                                const path_vec      &arg_pdb_filenames, ///< TODOCUMENT
                                                const bool          &arg_write_script,  ///< TODOCUMENT
                                                const bool          &arg_relabel_chain  ///< TODOCUMENT
                                                ) {
	// Read the PDBs
	pdb_list pdbs;
	pdbs.reserve(arg_pdb_filenames.size());
	for (const path &pdb_filename : arg_pdb_filenames) {
		pdb new_pdb;
		new_pdb.read_file(pdb_filename);
		pdbs.push_back(new_pdb);
	}

	write_superposed_pdb_to_file( arg_superposition, arg_filename, pdbs, arg_write_script, arg_relabel_chain );
}

/// \brief Save the specified superposition to the specified Boost Property Tree ptree
///
/// \relates superposition
void cath::sup::save_to_ptree(ptree               &arg_ptree,        ///< The ptree to which the superposition should be saved
                              const superposition &arg_superposition ///< The superposition to save to the ptree
                              ) {
	const auto transformations_key = string( "transformations" );
	arg_ptree.put_child( transformations_key, ptree{} );
	auto &transformations_ptree = arg_ptree.get_child( transformations_key );

	for (const size_t &index : irange( 0_z, arg_superposition.get_num_entries() ) ) {
		ptree transformation_ptree;
		transformation_ptree.put_child( "translation", make_ptree_of( arg_superposition.get_translation_of_index( index ) ) );
		transformation_ptree.put_child( "rotation",    make_ptree_of( arg_superposition.get_rotation_of_index   ( index ) ) );
		transformations_ptree.push_back( make_pair( "", transformation_ptree ) );
	}
}

/// \brief Make a new Boost Property Tree ptree representing the specified superposition
///
/// \relates superposition
ptree cath::sup::make_ptree_of(const superposition &arg_superposition ///< The superposition that the new ptree should represent
                               ) {
	ptree new_ptree;
	save_to_ptree( new_ptree, arg_superposition );
	return new_ptree;
}

/// \brief TODOCUMENT
///
/// \relates superposition
string cath::sup::to_json_string(const superposition &arg_superposition, ///< TODOCUMENT
                                 const bool          &arg_pretty_print   ///< TODOCUMENT
                                 ) {
	ostringstream json_ss;
	ptree temp_ptree;
	save_to_ptree( temp_ptree, arg_superposition );
	write_json( json_ss, temp_ptree, arg_pretty_print );
	return json_ss.str();
}
