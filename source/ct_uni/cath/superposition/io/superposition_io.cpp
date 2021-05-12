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

#include "superposition_io.hpp"

#include <filesystem>
#include <fstream>

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/range/adaptor/map.hpp>

#include <spdlog/spdlog.h>

#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/property_tree/make_ptree_of.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/outputter/superposition_output_options/superposition_output_options_block.hpp"
#include "cath/structure/structure_type_aliases.hpp"
#include "cath/superposition/superposition.hpp"

using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;
using namespace ::cath::opts;
using namespace ::cath::sup;
using namespace ::cath::sup::detail;

using ::boost::adaptors::map_values;
using ::boost::algorithm::any_of;
using ::boost::irange;
using ::boost::property_tree::ptree;
using ::std::filesystem::path;
using ::std::fixed;
using ::std::ofstream;
using ::std::ostream;
using ::std::setprecision;
using ::std::strerror;
using ::std::string;

const string superposition_io_consts::ENTRIES_KEY         = "entries";
const string superposition_io_consts::NAME_KEY            = "name";
const string superposition_io_consts::ROTATION_KEY        = "rotation";
const string superposition_io_consts::TRANSFORMATION_KEY  = "transformation";
const string superposition_io_consts::TRANSFORMATIONS_KEY = "transformations";
const string superposition_io_consts::TRANSLATION_KEY     = "translation";

/// \brief TODOCUMENT
///
/// \relates superposition
//double cath::get_rmsd_from_pairwise_superposition(const superposition &prm_superposition
//                                                  ) {
//	check_superposition_is_pairwise(prm_superposition);
//	return prm_superposition.get_rmsd_between_index_and_next(superposition::INDEX_OF_FIRST_IN_PAIRWISE_SUPERPOSITION);
//}

/// \brief Write a the superposition's essential data out to an XML superposition file
///
/// \relates superposition
///
/// \todo Improve the use of C++. Read:
///        - http://www.ibm.com/developerworks/xml/library/x-ctlbx/index.html
///        - http://www.grinninglizard.com/tinyxml2docs/index.html
///       (tinyxml is available as installable package on orengobuild64)
void cath::sup::write_xml_sup(ostream              &prm_ostream,       ///< TODOCUMENT
                              const superposition  &prm_superposition, ///< TODOCUMENT
                              const str_vec        &prm_ids            ///< TODOCUMENT
                              ) {
	if (prm_superposition.get_num_entries() != superposition::NUM_ENTRIES_IN_PAIRWISE_SUPERPOSITION) {
		BOOST_THROW_EXCEPTION(not_implemented_exception("Currently only able to output XML superposition for two entries"));
	}
	const size_t num_entries = prm_superposition.get_num_entries();
	if (prm_ids.size() != num_entries) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of IDs must match the number of entries in the superposition"));
	}

	prm_ostream << setprecision(6) << fixed;
	prm_ostream << "<?xml version=\"1.0\"?>\n";
	prm_ostream << "<root>\n";

	for (const size_t &entry_ctr : indices( num_entries ) ) {
		const coord centre_of_gravity = -prm_superposition.get_translation_of_index(entry_ctr);
		prm_ostream << "  <structure" << entry_ctr+1 << " id=\"" << prm_ids[entry_ctr];
		prm_ostream << "\">\n";
		prm_ostream << "    <centre x=\"" << centre_of_gravity.get_x();
		prm_ostream <<          "\" y=\"" << centre_of_gravity.get_y();
		prm_ostream <<          "\" z=\"" << centre_of_gravity.get_z();
		prm_ostream <<          "\" />\n";
		prm_ostream << "  </structure" << entry_ctr+1 << ">\n";
	}

	for (const size_t &rotation_ctr : irange( 1_z, num_entries ) ) {
		const rotation sup_rotn = prm_superposition.get_rotation_of_index(rotation_ctr);
		for (const size_t &dim_ctr_1 : indices( coord::NUM_DIMS ) ) {
			prm_ostream << "  <rotationmatrix name=\"row" << dim_ctr_1+1 << "\"";
			for (const size_t &dim_ctr_2 : indices( coord::NUM_DIMS ) ) {
				// NOTE: transposing the matrix to be consistent with legacy behaviour
				prm_ostream << " col" << dim_ctr_2+1 << "=\"" << sup_rotn.get_value(dim_ctr_2, dim_ctr_1) << "\"";
			}
			prm_ostream << " />\n";
		}
	}
	prm_ostream << "</root>\n";
}


/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::write_xml_sup_filename(const superposition  &prm_superposition, ///< TODOCUMENT
                                       const path           &prm_filename,      ///< TODOCUMENT
                                       const str_vec        &prm_ids            ///< TODOCUMENT
                                       ) {
	// A try block with the aim of catching any I/O errors
	try {
		ofstream xml_ostream = open_ofstream( prm_filename );
		write_xml_sup(xml_ostream, prm_superposition, prm_ids);
		xml_ostream.close();
	}
	catch (const boost::exception &ex) {
		throw;
	}
	catch (const std::exception & ex) {
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			  "Cannot write XML superposition file \""
			+ prm_filename.string()
			+ "\" ["
			+ ex.what()
			+ "] : "
			+ strerror( errno )
		));
	};
}

///// \brief TODOCUMENT
/////
///// \relates superposition
//ostream & cath::sup::write_superposed_pdb_to_ostream(ostream             &prm_os,            ///< TODOCUMENT
//                                                     const superposition &prm_superposition, ///< TODOCUMENT
//                                                     const pdb           &prm_pdb            ///< TODOCUMENT
//                                                     ) {
//	write_superposed_pdb_to_ostream(
//		prm_os,
//		prm_superposition,
//		prm_pdb,
//		0,
//		false
//	);
//	return prm_os;
//}

/// \brief TODOCUMENT
///
/// \relates superposition
ostream & cath::sup::write_superposed_pdb_to_ostream(ostream                    &prm_os,            ///< TODOCUMENT
                                                     const superposition        &prm_superposition, ///< TODOCUMENT
                                                     pdb                         prm_pdb,           ///< TODOCUMENT
                                                     const size_t               &prm_chain_index,   ///< TODOCUMENT
                                                     const chain_relabel_policy &prm_relabel_chain, ///< TODOCUMENT
                                                     const region_vec_opt       &prm_regions,       ///< Optional specification of regions to which the written records should be restricted
                                                     const pdb_write_mode       &prm_pdb_write_mode ///< Whether this is the only/last part of the PDB file
                                                     ) {
	// Translate PDB
	prm_pdb += prm_superposition.get_translation_of_index( prm_chain_index );

	// Rotate PDB
	prm_pdb.rotate( prm_superposition.get_rotation_of_index( prm_chain_index ) );

	if ( prm_relabel_chain == chain_relabel_policy::RELABEL ) {
		// Label structure with sensible chain label
		//
		// \todo Change this so that it does not relabel chains unless that is required
		//       (when writing a single old-style superposition PDB)
		//
		// *****
		// * THIS IS USED BY THE GENOME3D SUPERPOSITIONS SO ADD A TEST TO
		// * PRESERVE THAT BEHAVIOUR IF ALTERING THIS CODE
		// *****
		if ( prm_chain_index >= superposition::SUPERPOSITION_CHAIN_LABELS.size() ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Ran out of chain labels when relabelling chains"));
		}

		// If the PDB being written out has multiple chain_labels then setting the whole PDB's
		// chain_label to one value is probably a pretty bad idea so warn
		if ( has_multiple_chain_labels( prm_pdb ) ) {
			::spdlog::warn( "Overwriting chain labels whilst writing a PDB but this PDB contains multiple chain "
			                "labels, so this is probably a bad idea. Please consider different options for writing out "
			                "the superposed PDBs, such as --{}.",
			                superposition_output_options_block::PO_SUP_FILES_DIR );
		}

		prm_pdb.set_chain_label( superposition::SUPERPOSITION_CHAIN_LABELS[ prm_chain_index ] );
	}

	write_pdb_file( prm_os, prm_pdb, prm_regions, prm_pdb_write_mode );

	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates superposition
ostream & cath::sup::write_superposed_pdbs_to_ostream(ostream                      &prm_os,             ///< TODOCUMENT
                                                      const superposition          &prm_superposition,  ///< TODOCUMENT
                                                      const pdb_list                prm_pdbs,           ///< TODOCUMENT
                                                      const sup_pdbs_script_policy &prm_script_policy,  ///< TODOCUMENT
                                                      const chain_relabel_policy   &prm_relabel_chain,  ///< TODOCUMENT
                                                      const region_vec_opt         &prm_regions         ///< Optional specification of regions to which the written records should be restricted
                                                      ) {
	// Sanity check inputs
	const size_t num_entries = prm_superposition.get_num_entries();
	if (num_entries != prm_pdbs.size()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot superpose PDBs because the number of PDBs does not match the number of entries in the superposition"));
	}
	if (num_entries > superposition::SUPERPOSITION_CHAIN_LABELS.size()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot superpose PDBs because there aren't enough chain labels"));
	}

	if ( prm_script_policy == sup_pdbs_script_policy::WRITE_RASMOL_SCRIPT ) {
		prm_os << R"(#!rasmol -script
zap  
load inline
wireframe off
select all
cartoon
select *
color chain
exit
)";
	}

	// Translate each PDB to midpoint, based on CoG of equivalent positions
	// (This is important for structures with low overlap)
	for (const size_t pdb_ctr : indices( num_entries ) ) {
		write_superposed_pdb_to_ostream(
			prm_os,
			prm_superposition,
			prm_pdbs[pdb_ctr],
			pdb_ctr,
			prm_relabel_chain,
			prm_regions,
			( ( pdb_ctr + 1 == num_entries ) ? pdb_write_mode::ONLY_OR_LAST_PDB : pdb_write_mode::MORE_TO_FOLLOW )
		);
	}

	return prm_os;
}

///// TODOCUMENT
/////
///// \relates superposition
/////
///// \todo Refactor out common file-opening/exception-handling functionality from these file readers/writers
//void cath::sup::write_superposed_pdb_to_file(const superposition &prm_superposition, ///< TODOCUMENT
//                                             const path          &prm_filename,      ///< TODOCUMENT
//                                             const pdb           &prm_pdb            ///< TODOCUMENT
//                                             ) {
//	ofstream superposition_ostream;
//	open_ofstream(superposition_ostream, prm_filename);
//
//	// Try here to catch any I/O exceptions
//	try {
//		write_superposed_pdb_to_ostream(
//			superposition_ostream,
//			prm_superposition,
//			prm_pdb,
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
//		BOOST_THROW_EXCEPTION(runtime_error_exception(
//			  "Cannot write superposition to file \""
//			+ prm_filename.string()
//			+ "\" ["
//			+ ex.what()
//			+ "] : "
//			+ strerror( errno )
//		));
//	};
//}

/// TODOCUMENT
///
/// \relates superposition
///
/// \todo Refactor out common file-opening/exception-handling functionality from these file readers/writers
void cath::sup::write_superposed_pdb_to_file(const superposition        &prm_superposition,  ///< TODOCUMENT
                                             const path                 &prm_filename,       ///< TODOCUMENT
                                             const pdb                  &prm_pdb,            ///< TODOCUMENT
                                             const size_t               &prm_chain_index,    ///< TODOCUMENT
                                             const chain_relabel_policy &prm_relabel_chain,  ///< TODOCUMENT
                                             const region_vec_opt       &prm_regions         ///< Optional specification of regions to which the written records should be restricted
                                             ) {
	ofstream superposition_ostream = open_ofstream( prm_filename);

	// Try here to catch any I/O exceptions
	try {
		write_superposed_pdb_to_ostream(
			superposition_ostream,
			prm_superposition,
			prm_pdb,
			prm_chain_index,
			prm_relabel_chain,
			prm_regions
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
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			  "Cannot write superposition to file \""
			+ prm_filename.string()
			+ "\" ["
			+ ex.what()
			+ "] : "
			+ strerror( errno )
		));
	};
}

/// TODOCUMENT
///
/// \relates superposition
///
/// \todo Refactor out common file-opening/exception-handling functionality from these file readers/writers
void cath::sup::write_superposed_pdb_to_file(const superposition          &prm_superposition,  ///< TODOCUMENT
                                             const path                   &prm_filename,       ///< TODOCUMENT
                                             const pdb_list               &prm_pdbs,           ///< TODOCUMENT
                                             const sup_pdbs_script_policy &prm_script_policy,  ///< TODOCUMENT
                                             const chain_relabel_policy   &prm_relabel_chain,  ///< TODOCUMENT
                                             const region_vec_opt         &prm_regions         ///< Optional specification of regions to which the written records should be restricted
                                             ) {
	ofstream superposition_ostream = open_ofstream( prm_filename );

	// Try here to catch any I/O exceptions
	try {
		write_superposed_pdbs_to_ostream(
			superposition_ostream,
			prm_superposition,
			prm_pdbs,
			prm_script_policy,
			prm_relabel_chain,
			prm_regions
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
		BOOST_THROW_EXCEPTION(runtime_error_exception(
			  "Cannot write superposition to file \""
			+ prm_filename.string()
			+ "\" ["
			+ ex.what()
			+ "] : "
			+ strerror( errno )
		));
	};
}

/// \brief TODOCUMENT
///
/// \relates superposition
void cath::sup::write_superposed_pdb_from_files(const superposition          &prm_superposition,  ///< TODOCUMENT
                                                const path                   &prm_filename,       ///< TODOCUMENT
                                                const path_vec               &prm_pdb_filenames,  ///< TODOCUMENT
                                                const sup_pdbs_script_policy &prm_script_policy,  ///< TODOCUMENT
                                                const chain_relabel_policy   &prm_relabel_chain,  ///< TODOCUMENT
                                                const region_vec_opt         &prm_regions         ///< Optional specification of regions to which the written records should be restricted
                                                ) {
	// Read the PDBs
	pdb_list pdbs;
	pdbs.reserve(prm_pdb_filenames.size());
	for (const path &pdb_filename : prm_pdb_filenames) {
		pdb new_pdb;
		new_pdb.read_file(pdb_filename);
		pdbs.push_back(new_pdb);
	}

	write_superposed_pdb_to_file(
		prm_superposition,
		prm_filename,
		pdbs,
		prm_script_policy,
		prm_relabel_chain,
		prm_regions
	);
}

/// \brief Build a superposition from a superposition-populated ptree
///
/// \relates superposition
superposition cath::sup::superposition_from_ptree(const ptree &prm_ptree ///< The ptree from which the superposition should be read
                                                  ) {
	// Sanity check the ptree
	const auto tranformations      = prm_ptree.get_child( superposition_io_consts::TRANSFORMATIONS_KEY );
	const auto num_transformations = tranformations.size();
	if ( num_transformations != tranformations.count( "" ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse superposition from property_tree with any unrecognised transformations fields"));
	}
	const auto transformation_entry_is_invalid = [] (const ptree &x) {
		return (
			   x.size ()                                           != 2
			|| x.count( superposition_io_consts::TRANSLATION_KEY ) != 1
			|| x.count( superposition_io_consts::ROTATION_KEY    ) != 1
		);
	};
	if ( any_of( tranformations | map_values, transformation_entry_is_invalid ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse superposition from property_tree with invalid transformation entry"));
	}

	// Read the translations
	const auto translations = transform_build<coord_vec>(
		tranformations | map_values,
		[] (const ptree &x) {
			return coord_from_ptree( x.get_child( superposition_io_consts::TRANSLATION_KEY ) );
		}
	);

	// Parse the rotations
	const auto rotations = transform_build<rotation_vec>(
		tranformations | map_values,
		[] (const ptree &x) {
			return rotation_from_ptree( x.get_child( superposition_io_consts::ROTATION_KEY ) );
		}
	);

	// Return a superposition of these translations and rotations
	return { translations, rotations };
}

/// \brief Save the specified superposition to the specified Boost Property Tree ptree
///
/// \relates superposition
void cath::sup::save_to_ptree(ptree               &prm_ptree,        ///< The ptree to which the superposition should be saved
                              const superposition &prm_superposition ///< The superposition to save to the ptree
                              ) {
	const auto transformations_key = string( superposition_io_consts::TRANSFORMATIONS_KEY );
	prm_ptree.put_child( transformations_key, ptree{} );
	auto &transformations_ptree = prm_ptree.get_child( transformations_key );

	for (const size_t &index : indices( prm_superposition.get_num_entries() ) ) {
		ptree transformation_ptree;
		transformation_ptree.put_child( superposition_io_consts::TRANSLATION_KEY, make_ptree_of( prm_superposition.get_translation_of_index( index ) ) );
		transformation_ptree.put_child( superposition_io_consts::ROTATION_KEY,    make_ptree_of( prm_superposition.get_rotation_of_index   ( index ) ) );
		transformations_ptree.push_back( make_pair( "", transformation_ptree ) );
	}
}
