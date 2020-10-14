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

#include "superposition_context.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include "alignment/alignment_context.hpp"
#include "chopping/region/region.hpp"
#include "common/algorithm/append.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "common/property_tree/make_ptree_of.hpp"
#include "file/options/data_dirs_options_block.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "file/pdb/proximity_calculator.hpp"
#include "structure/structure_type_aliases.hpp"
#include "superposition/io/superposition_io.hpp"
#include "superposition/superposition_content_spec.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;
using namespace ::cath::opts;
using namespace ::cath::sup;
using namespace ::cath::sup::detail;

using ::boost::adaptors::map_values;
using ::boost::adaptors::transformed;
using ::boost::algorithm::any_of;
using ::boost::log::trivial::severity_level;
using ::boost::make_optional;
using ::boost::none;
using ::boost::property_tree::ptree;
using ::boost::range::for_each;
using ::std::pair;
using ::std::string;

/// \brief Ctor for superposition_context
superposition_context::superposition_context(superposition  prm_superposition, ///< TODOCUMENT
                                             strucs_context prm_context        ///< TODOCUMENT
                                             ) : the_superposition { std::move( prm_superposition ) },
                                                 context           { std::move( prm_context       ) } {
}

/// \brief Ctor for superposition_context
superposition_context::superposition_context(superposition  prm_superposition, ///< TODOCUMENT
                                             strucs_context prm_context,       ///< TODOCUMENT
                                             alignment      prm_alignment      ///< TODOCUMENT
                                             ) : the_superposition { std::move( prm_superposition ) },
                                                 context           { std::move( prm_context       ) },
                                                 any_alignment     { std::move( prm_alignment     ) } {
}

/// \brief Ctor for superposition_context
superposition_context::superposition_context(superposition      prm_superposition, ///< TODOCUMENT
                                             pdb_list           prm_pdbs,          ///< TODOCUMENT
                                             name_set_list      prm_name_sets,     ///< TODOCUMENT
                                             region_vec_opt_vec prm_regions        ///< The key regions of the structures
                                             ) : the_superposition { std::move( prm_superposition )   },
                                                 context           {
                                                 	std::move( prm_pdbs      ),
                                                 	std::move( prm_name_sets ),
                                                 	std::move( prm_regions   )
                                                 } {
}

/// \brief Ctor for superposition_context
superposition_context::superposition_context(superposition      prm_superposition, ///< TODOCUMENT
                                             pdb_list           prm_pdbs,          ///< TODOCUMENT
                                             name_set_list      prm_name_sets,     ///< TODOCUMENT
                                             region_vec_opt_vec prm_regions,       ///< The key regions of the structures
                                             alignment          prm_alignment      ///< TODOCUMENT
                                             ) : the_superposition { std::move( prm_superposition )   },
                                                 context           {
                                                 	std::move( prm_pdbs      ),
                                                 	std::move( prm_name_sets ),
                                                 	std::move( prm_regions   )
                                                 },
                                                 any_alignment     { std::move( prm_alignment     )   } {
}

/// \brief TODOCUMENT
const superposition & superposition_context::get_superposition() const {
	return the_superposition;
}

/// \brief TODOCUMENT
const strucs_context & superposition_context::get_strucs_context() const {
	return context;
}

/// \brief TODOCUMENT
bool superposition_context::has_alignment() const {
	return static_cast<bool>( any_alignment );
}

/// \brief TODOCUMENT
const alignment & superposition_context::get_alignment() const {
	if ( ! has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to get alignment from superposition_context that doesn't contain one"));
	}
	return *any_alignment;
}

/// \brief Setter for the PDBs
///
/// \pre prm_pdbs.size() == get_num_entries( *this ) else this throws an invalid_argument_exception
superposition_context & superposition_context::set_pdbs(const pdb_list &prm_pdbs ///< The PDBs to set
                                                        ) {
	if ( prm_pdbs.size() != get_num_entries( *this ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Unable to load "                                           + ::std::to_string( prm_pdbs.size()          )
			+ " pdbs into superposition context of superposition with " + ::std::to_string( get_num_entries( *this ) )
			+ " entries"
		));
	}
	context.set_pdbs( prm_pdbs );
	return *this;
}

/// \brief TODOCUMENT
///
/// \relates superposition_context
const pdb_list & cath::sup::get_pdbs(const superposition_context &prm_supn_context ///< TODOCUMENT
                                     ) {
	return prm_supn_context.get_strucs_context().get_pdbs();
}

/// \brief TODOCUMENT
///
/// \relates superposition_context
const name_set_list & cath::sup::get_name_sets(const superposition_context &prm_supn_context ///< TODOCUMENT
                                               ) {
	return prm_supn_context.get_strucs_context().get_name_sets();
}

/// \brief Getter for the specification of the key regions of the structures
///
/// \relates superposition_context
const region_vec_opt_vec & cath::sup::get_regions(const superposition_context &prm_supn_context ///< TODOCUMENT
                                                  ) {
	return prm_supn_context.get_strucs_context().get_regions();
}

///// \brief TODOCUMENT
/////
///// \relates superposition_context
//bool cath::sup::operator==(const superposition_context &prm_sup_con_a, ///< TODOCUMENT
//                           const superposition_context &prm_sup_con_b  ///< TODOCUMENT
//                           ) {
//}

///// \brief TODOCUMENT
/////
///// \relates superposition_context
//ostream & cath::sup::operator<<(ostream                     &prm_os,     ///< TODOCUMENT
//                                const superposition_context &prm_sup_con ///< TODOCUMENT
//                                ) {
//}

/// \brief Get a copy of the PDBs in the specified superposition_context, restricted to its regions
///
/// \relates superposition_context
pdb_list cath::sup::get_restricted_pdbs(const superposition_context &prm_superposition_context ///< The superposition_context to query
                                        ) {
	return get_restricted_pdbs( prm_superposition_context.get_strucs_context() );
}

/// \brief Get a PDB with the appropriate parts of the specified PDB according to the specified
///        regions (or none for all) and specified superposition_content_spec
///
/// At the moment, if the supn_regions_context is IN_CHAIN or IN_PDB, this will still only
/// ligand atoms based on the specified region(s).
/// This could be altered (or made configurable in the future).
///
/// If the regions are not specified (ie none), the whole PDB is returned.
///
/// This does a similar job to get_regions_limited_pdb() but adds in the extra
/// handling of the superposition_content_spec.
pdb cath::sup::get_supn_content_pdb(const pdb                        &prm_pdb,         ///< The source PDB
                                    const region_vec_opt             &prm_regions,     ///< The key regions of the structure
                                    const superposition_content_spec &prm_content_spec ///< The specification of what should be included in the superposition
                                    ) {
	if ( ! prm_regions ) {
		return prm_pdb;
	}

	// From http://www.proteopedia.org/wiki/index.php/Hydrogen_bonds:
	//
	// > Jeffrey[1] (page 12) categorizes hbonds with donor-acceptor distances of
	// > 2.2-2.5 Å as "strong, mostly covalent",
	// > 2.5-3.2 Å as "moderate, mostly electrostatic", and
	// > 3.2-4.0 Å as "weak, electrostatic".
	//
	// ...where reference [1] is:
	// Jeffrey, George A., An introduction to hydrogen bonding, Oxford University Press, 1997.
	constexpr double EXTEND_DIST = 4.2;
	const doub_opt &include_dna_within_distance     = prm_content_spec.get_include_dna_within_distance();
	const doub_opt &include_organic_within_distance = prm_content_spec.get_include_organic_within_distance();

	const auto      expanded_regions = get_regions_expanded_for_context( *prm_regions, prm_content_spec );

	const proximity_calculator prox_calc{ prm_pdb, *prm_regions };

	const chain_label_set nearby_dna_chains = nearby_dna_rna_chain_labels(
		prm_pdb,
		prox_calc,
		include_dna_within_distance
	);
	const pdb_residue_vec nearby_post_ter_coords = get_nearby_post_ter_res_atoms(
		prm_pdb,
		prox_calc,
		include_organic_within_distance,
		EXTEND_DIST
	);

	// Form a list of regions that includes the originals,
	// plus any DNA chains that have been found nearby
	const region_vec_opt all_regions = expanded_regions
		? ::boost::make_optional( append_copy(
			*expanded_regions,
			transform_build<region_vec>(
				nearby_dna_chains,
				[] (const chain_label &x) { return region{ x }; }
			)
		) )
		: none;

	// Build a PDB from the information and return it
	pdb result_pdb = get_regions_limited_pdb(
		all_regions,
		prm_pdb
	);
	result_pdb.set_post_ter_residues( nearby_post_ter_coords );
	return result_pdb;
}

/// \brief Get PDBs with the appropriate parts of the specified superposition_context's PDBs according to its
///        regions and the specified superposition_content_spec
///
/// \relates superposition_context
pdb_list cath::sup::get_supn_content_pdbs(const superposition_context      &prm_superposition_context, ///< The superposition_context containing the PDBs and regions to extract
                                          const superposition_content_spec &prm_content_spec           ///< The specification of what should be included in the superposition 
                                          ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return pdb_list{
		transform_build<pdb_vec>(
			get_pdbs   ( prm_superposition_context ),
			get_regions( prm_superposition_context ),
			[&] (const pdb &the_pdb, const region_vec_opt &the_regions) {
				return get_supn_content_pdb( the_pdb, the_regions, prm_content_spec );
			}
		)
	};
}

/// \brief Return a copy of the specified superposition_context in which the specified PDBs have been set
///
/// \relates superposition_context
superposition_context cath::sup::set_pdbs_copy(superposition_context  prm_sup_con, ///< The superposition_context from which a copy should be taken, altered and returned
                                               const pdb_list        &prm_pdbs     ///< The PDBs to set on the copy of the superposition_context
                                               ) {
	prm_sup_con.set_pdbs( prm_pdbs );
	return prm_sup_con;
}

/// \brief Get the number of entries in the specified superposition_context
///
/// \relates superposition_context
size_t cath::sup::get_num_entries(const superposition_context &prm_superposition_context ///< The superposition_context to query
                                  ) {
	return prm_superposition_context.get_superposition().get_num_entries();
}

/// \brief Load the specified superposition_context's PDBs using its names and the specified data_dirs_options_block
///
/// \relates superposition_context
void cath::sup::load_pdbs_from_names(superposition_context &prm_supn_context, ///< The superposition_context for which the PDBs should be loaded based on its names
                                     const data_dirs_spec  &prm_data_dirs     ///< The data_dirs_options_block with which to convert names into PDB filenames
                                     ) {
	prm_supn_context.set_pdbs(
		pdb_list{ transform_build<pdb_vec>(
			get_name_sets( prm_supn_context ),
			[&] (const name_set &the_name_set) {
				return read_pdb_file( find_file( prm_data_dirs, data_file::PDB, the_name_set ) );
			}
		) }
	);
}

/// \brief Return a copy of the specified superposition_context with the PDBs loaded using its names and the specified data_dirs_options_block
///
/// \relates superposition_context
superposition_context cath::sup::load_pdbs_from_names_copy(superposition_context  prm_supn_context, ///< The superposition_context from which the copy should be taken with the PDBs loaded based on its names
                                                           const data_dirs_spec  &prm_data_dirs     ///< The data_dirs_options_block with which to convert names into PDB filenames
                                                           ) {
	load_pdbs_from_names( prm_supn_context, prm_data_dirs );
	return prm_supn_context;
}

/// \brief TODOCUMENT
///
/// \relates superposition_context
///
/// \relatesalso alignment_context
alignment_context cath::sup::make_restricted_alignment_context(const superposition_context &prm_superposition_context ///< TODOCUMENT
                                                               ) {
	if ( ! prm_superposition_context.has_alignment() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("WARNING: Unable to extract alignment from superposition_context with no alignment"));
	}
	return {
		prm_superposition_context.get_alignment(),
		get_restricted_pdbs( prm_superposition_context ),
		get_name_sets      ( prm_superposition_context ),
		get_regions        ( prm_superposition_context )
	};
}

/// \brief  Build a coord from a superposition_context-populated ptree
///
/// \relates superposition_context
superposition_context cath::sup::superposition_context_from_ptree(const ptree &prm_ptree ///< The ptree from which the superposition_context should be read
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
	if ( prm_ptree.size() != 1 || prm_ptree.count( superposition_io_consts::ENTRIES_KEY ) != 1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse a superposition_context from ptree data that doesn't have one entries key and no other keys"));
	}
	const auto entries = prm_ptree.get_child( superposition_io_consts::ENTRIES_KEY );

	// Sanity check the ptree [ Step 2: check that all entries have empty keys ]
	if ( entries.size() != entries.count( "" ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot parse a superposition_context from ptree data whose entries have non-empty"));
	}

	// Sanity check the ptree [ Step 3: check that all values contain exactly two keys, name and transformation ]
	if ( any_of( entries | map_values, entry_is_invalid ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Cannot parse a superposition_context from ptree data whose entries don't contain exactly two entries: name and transformation"));
	}

	// Read the names
	const name_set_list names { transform_build<name_set_vec>(
		entries | map_values,
		[] (const ptree &x) {
			return name_set{ x.get<string>( superposition_io_consts::NAME_KEY ) };
		}
	) };

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
		superposition{ translations, rotations },
		pdb_list{ pdb_vec{ names.size() } },
		names,
		region_vec_opt_vec{ names.size() }
	};
}

/// \brief TODOCUMENT
///
/// At present, this stores the names and the superposition but does nothing
/// with the alignment or the PDBs
///
/// \relates superposition_context
void cath::sup::save_to_ptree(ptree                       &prm_ptree,      ///< TODOCUMENT
                              const superposition_context &prm_sup_context ///< TODOCUMENT
                              ) {
	if ( prm_sup_context.has_alignment() ) {
		BOOST_LOG_TRIVIAL( warning ) << "Whilst converting a superposition_context to JSON, its alignment will be ignored because that is not currently supported";
	}

	const auto supn_ptree          = make_ptree_of( prm_sup_context.get_superposition() );
	const auto trans_ptrees        = supn_ptree.get_child( superposition_io_consts::TRANSFORMATIONS_KEY );

	prm_ptree.put_child( superposition_io_consts::ENTRIES_KEY, ptree{} );
	auto &entries_ptree = prm_ptree.get_child( superposition_io_consts::ENTRIES_KEY );

	for_each(
		get_supn_json_names( get_name_sets( prm_sup_context ) ),
		trans_ptrees,
		[&] (const string &name, const pair<string, ptree> &trans_ptree) {
			ptree entry_ptree;
			entry_ptree.put      ( superposition_io_consts::NAME_KEY,           name               );
			entry_ptree.put_child( superposition_io_consts::TRANSFORMATION_KEY, trans_ptree.second );
			entries_ptree.push_back( make_pair( "", entry_ptree ) );
		}
	);
}
