/// \file
/// \brief The superposition_context class definitions

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

#include "superposition_context.h"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include "alignment/alignment_context.h"
#include "common/algorithm/transform_build.h"
#include "exception/invalid_argument_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "options/options_block/data_dirs_options_block.h"
#include "structure/structure_type_aliases.h"
#include "superposition/io/superposition_io.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::opts;
using namespace cath::sup;
using namespace cath::sup::detail;
using namespace std;

using boost::adaptors::map_values;
using boost::adaptors::transformed;
using boost::algorithm::any_of;
using boost::filesystem::path;
using boost::log::trivial::severity_level;
using boost::property_tree::json_parser::write_json;
using boost::property_tree::ptree;
using boost::range::for_each;

/// \brief TODOCUMENT
superposition_context::superposition_context(const pdb_list      &arg_pdbs,         ///< TODOCUMENT
                                             const str_vec       &arg_names,        ///< TODOCUMENT
                                             const superposition &arg_superposition ///< TODOCUMENT
                                             ) : pdbs(arg_pdbs),
                                                 names(arg_names),
                                                 the_superposition(arg_superposition) {
}


/// \brief TODOCUMENT
superposition_context::superposition_context(const pdb_list      &arg_pdbs,          ///< TODOCUMENT
                                             const str_vec       &arg_names,         ///< TODOCUMENT
                                             const superposition &arg_superposition, ///< TODOCUMENT
                                             const alignment     &arg_alignment      ///< TODOCUMENT
                                             ) : pdbs(arg_pdbs),
                                                 names(arg_names),
                                                 the_superposition(arg_superposition),
                                                 any_alignment(1, arg_alignment) {

}

/// \brief TODOCUMENT
const pdb_list & superposition_context::get_pdbs_cref() const {
	return pdbs;
}

/// \brief TODOCUMENT
const str_vec & superposition_context::get_names_cref() const {
	return names;
}

/// \brief TODOCUMENT
const superposition & superposition_context::get_superposition_cref() const {
	return the_superposition;
}

/// \brief TODOCUMENT
bool superposition_context::has_alignment() const {
	return static_cast<bool>( any_alignment );
}

/// \brief TODOCUMENT
const alignment & superposition_context::get_alignment_cref() const {
	if ( ! has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get alignment from superposition_context that doesn't contain one"));
	}
	return *any_alignment;
}

/// \brief Setter for the PDBs
///
/// \pre arg_pdbs.size() == get_num_entries( *this ) else this throws an invalid_argument_exception
void superposition_context::set_pdbs(const pdb_list &arg_pdbs ///< The PDBs to set
                                     ) {
	if ( arg_pdbs.size() != get_num_entries( *this ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to load "                                           + to_string( arg_pdbs.size()          )
			+ " pdbs into superposition context of superposition with " + to_string( get_num_entries( *this ) )
			+ " entries"
		));
	}
	pdbs = arg_pdbs;
}

///// \brief TODOCUMENT
/////
///// \relates superposition_context
//bool cath::sup::operator==(const superposition_context &arg_sup_con_a, ///< TODOCUMENT
//                           const superposition_context &arg_sup_con_b  ///< TODOCUMENT
//                           ) {
//}

///// \brief TODOCUMENT
/////
///// \relates superposition_context
//ostream & cath::sup::operator<<(ostream                     &arg_os,     ///< TODOCUMENT
//                                const superposition_context &arg_sup_con ///< TODOCUMENT
//                                ) {
//}

/// \brief Get the number of entries in the specified superposition_context
///
/// \relates superposition_context
size_t cath::sup::get_num_entries(const superposition_context &arg_superposition_context ///< The superposition_context to query
                                  ) {
	return arg_superposition_context.get_superposition_cref().get_num_entries();
}

/// \brief Load the specified superposition_context's PDBs using its names and the specified data_dirs_options_block
///
/// \relates superposition_context
void cath::sup::load_pdbs_from_names(superposition_context         &arg_supn_context, ///< The superposition_context for which the PDBs should be loaded based on its names
                                     const data_dirs_options_block &arg_data_dirs     ///< The data_dirs_options_block with which to convert names into PDB filenames
                                     ) {
	arg_supn_context.set_pdbs(
		transform_build<pdb_vec>(
			arg_supn_context.get_names_cref(),
			[&] (const string &name) {
				return read_pdb_file( find_file( arg_data_dirs, data_file::PDB, name ) );
			}
		)
	);
}

/// \brief Return a copy of the specified superposition_context with the PDBs loaded using its names and the specified data_dirs_options_block
///
/// \relates superposition_context
superposition_context cath::sup::load_pdbs_from_names_copy(superposition_context          arg_supn_context, ///< The superposition_context from which the copy should be taken with the PDBs loaded based on its names
                                                           const data_dirs_options_block &arg_data_dirs     ///< The data_dirs_options_block with which to convert names into PDB filenames
                                                           ) {
	load_pdbs_from_names( arg_supn_context, arg_data_dirs );
	return arg_supn_context;
}

/// \brief TODOCUMENT
///
/// \relates superposition_context
///
/// \relatesalso alignment_context
alignment_context cath::sup::make_alignment_context(const superposition_context &arg_superposition_context ///< TODOCUMENT
                                                    ) {
	if ( ! arg_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("WARNING: Unable to extract alignment from superposition_context with no alignment"));
	}
	return alignment_context(
		arg_superposition_context.get_pdbs_cref(),
		arg_superposition_context.get_names_cref(),
		arg_superposition_context.get_alignment_cref()
	);
}

/// \brief  Build a coord from a superposition_context-populated ptree
///
/// \relates superposition_context
superposition_context cath::sup::superposition_context_from_ptree(const ptree &arg_ptree ///< The ptree from which the superposition_context should be read
                                                                  ) {
	// Define a lambda for checking whether an entry ptree is invalid
	const auto entry_is_invalid = [] (const ptree &x) {
		return ( x.size() != 2
			  || x.count( superposition_io_consts::NAME_KEY           ) != 1
			  || x.count( superposition_io_consts::TRANSFORMATION_KEY ) != 1 );
	};

	// Define a lambda for returning an entry ptree's transformation ptree child
	const auto get_transformation_child = [] (const ptree &x) {
		return x.get_child( superposition_io_consts::TRANSFORMATION_KEY );
	};

	// Sanity check the ptree [ Step 1: check there's one key, which is entries ]
	if ( arg_ptree.size() != 1 || arg_ptree.count( superposition_io_consts::ENTRIES_KEY ) != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse a superposition_context from ptree data that doesn't have one entries key and no other keys"));
	}
	const auto entries = arg_ptree.get_child( superposition_io_consts::ENTRIES_KEY );

	// Sanity check the ptree [ Step 2: check that all entries have empty keys ]
	if ( entries.size() != entries.count( "" ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse a superposition_context from ptree data whose entries have non-empty"));
	}

	// Sanity check the ptree [ Step 3: check that all values contain exactly two keys, name and transformation ]
	if ( any_of( entries | map_values, entry_is_invalid ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Cannot parse a superposition_context from ptree data whose entries don't contain exactly two entries: name and transformation"));
	}

	// Read the names
	const auto names = transform_build<str_vec>(
		entries | map_values,
		[] (const ptree &x) {
			return x.get<string>( superposition_io_consts::NAME_KEY );
		}
	);

	// Read the translations
	const auto translations = transform_build<coord_vec>(
		entries | map_values | transformed( get_transformation_child ),
		[] (const ptree &x) {
			return coord_from_ptree( x.get_child( superposition_io_consts::TRANSLATION_KEY ) );
		}
	);

	// Parse the rotations
	const auto rotations = transform_build<rotation_vec>(
		entries | map_values | transformed( get_transformation_child ),
		[] (const ptree &x) {
			return rotation_from_ptree( x.get_child( superposition_io_consts::ROTATION_KEY ) );
		}
	);

	// Return a superposition_context built from the parsed data
	return {
		pdb_list{ pdb_vec{ names.size() } },
		names,
		superposition{ translations, rotations }
	};
}

/// \brief TODOCUMENT
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
void cath::sup::save_to_ptree(ptree                       &arg_ptree,      ///< TODOCUMENT
                              const superposition_context &arg_sup_context ///< TODOCUMENT
                              ) {
	if ( arg_sup_context.has_alignment() ) {
		BOOST_LOG_TRIVIAL( warning ) << "Whilst converting a superposition_context to JSON, its alignment will be ignored because that is not currently supported";
	}

	const auto supn_ptree          = make_ptree_of( arg_sup_context.get_superposition_cref() );
	const auto trans_ptrees        = supn_ptree.get_child( superposition_io_consts::TRANSFORMATIONS_KEY );

	arg_ptree.put_child( superposition_io_consts::ENTRIES_KEY, ptree{} );
	auto &entries_ptree = arg_ptree.get_child( superposition_io_consts::ENTRIES_KEY );

	for_each(
		arg_sup_context.get_names_cref(),
		trans_ptrees,
		[&] (const string &name, const pair<string, ptree> &trans_ptree) {
			ptree entry_ptree;
			entry_ptree.put      ( superposition_io_consts::NAME_KEY,           name               );
			entry_ptree.put_child( superposition_io_consts::TRANSFORMATION_KEY, trans_ptree.second );
			entries_ptree.push_back( make_pair( "", entry_ptree ) );
		}
	);
}

/// \brief TODOCUMENT
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
ptree cath::sup::make_ptree_of(const superposition_context &arg_sup_context ///< TODOCUMENT
                               ) {
	ptree new_ptree;
	save_to_ptree( new_ptree, arg_sup_context );
	return new_ptree;
}

/// \brief Build a superposition_context from a JSON string (via a ptree)
///
/// \relates superposition_context
superposition_context cath::sup::superposition_context_from_json_string(const string &arg_json_string ///< The JSON string from which the superposition_context should be read
                                                                        ) {
	ptree tree;
	istringstream in_ss( arg_json_string );
	read_json( in_ss, tree);
	return superposition_context_from_ptree( tree );
}

/// \brief Create a JSON string to represent the specified superposition
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
string cath::sup::to_json_string(const superposition_context &arg_sup_context, ///< The superposition_context to represent in the JSON string
                                 const bool                  &arg_pretty_print ///< Whether to use whitespace (including line breaks) in the JSON to make it more human-readable
                                 ) {
	ostringstream json_ss;
	ptree temp_ptree;
	save_to_ptree( temp_ptree, arg_sup_context );
	write_json( json_ss, temp_ptree, arg_pretty_print );
	return json_ss.str();
}
