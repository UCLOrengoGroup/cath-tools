/// \file
/// \brief The pdb_list class definitions

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

#include "pdb_list.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/range/combine.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "file/name_set/name_set_list.hpp"
#include "file/pdb/backbone_complete_indices.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "file/pdb/protein_info.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::filesystem::path;
using boost::range::combine;

/// \brief Ctor from a vector<pdb>
pdb_list::pdb_list(pdb_vec arg_pdbs ///< The pdbs from which this pdb_list should be constructed
                   ) : pdbs{ std::move( arg_pdbs ) } {
}

/// \brief TODOCUMENT
void pdb_list::push_back(const pdb &arg_pdb ///< TODOCUMENT
                         ) {
	pdbs.push_back(arg_pdb);
}

/// \brief TODOCUMENT
void pdb_list::reserve(const size_t &arg_size ///< TODOCUMENT
                       ) {
	pdbs.reserve(arg_size);
}

/// \brief TODOCUMENT
size_t pdb_list::size() const {
	return pdbs.size();
}

/// \brief TODOCUMENT
bool pdb_list::empty() const {
	return pdbs.empty();
}

/// \brief TODOCUMENT
pdb & pdb_list::operator[](const size_t &arg_index ///< TODOCUMENT
                           ) {
	return pdbs[arg_index];
}

/// \brief TODOCUMENT
const pdb & pdb_list::operator[](const size_t &arg_index ///< TODOCUMENT
                                 ) const {
	return pdbs[arg_index];
}

///// \brief TODOCUMENT
//pdb_list::iterator pdb_list::begin() {
//	return std::begin( pdbs );
//}
///// \brief TODOCUMENT
//pdb_list::iterator pdb_list::end() {
//	return std::end( pdbs );
//}
/// \brief TODOCUMENT
pdb_list::const_iterator pdb_list::begin() const {
	return common::cbegin( pdbs );
}
/// \brief TODOCUMENT
pdb_list::const_iterator pdb_list::end() const {
	return common::cend( pdbs );
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
pdb_list cath::file::read_pdb_files(const path_vec &arg_paths ///< TODOCUMENT
                                    ) {
	return make_pdb_list(
		transform_build<pdb_vec>(
			arg_paths,
			[] (const path &x) { return read_pdb_file( x ); }
		)
	);
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
pdb_list cath::file::make_pdb_list(const pdb_vec &arg_pdbs ///< TODOCUMENT
                                   ) {
	pdb_list new_pdb_list;
	new_pdb_list.reserve( arg_pdbs.size() );
	for (const pdb &the_pdb : arg_pdbs) {
		new_pdb_list.push_back( the_pdb );
	}
	return new_pdb_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
pdb_list cath::file::pdb_list_of_backbone_complete_subset_pdbs(const pdb_list        &arg_pdb_list, ///< TODOCUMENT
                                                               const ostream_ref_opt &arg_ostream   ///< An optional reference to an ostream to which any logging should be sent
                                                               ) {
	pdb_list new_pdb_list;
	new_pdb_list.reserve( arg_pdb_list.size() );
	for (const pdb &the_pdb : arg_pdb_list) {
		new_pdb_list.push_back( backbone_complete_subset_of_pdb( the_pdb, arg_ostream ).first );
	}
	return new_pdb_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
pdb_list cath::file::pdb_list_of_backbone_complete_region_limited_subset_pdbs(const pdb_list           &arg_pdb_list, ///< TODOCUMENT
                                                                              const region_vec_opt_vec &arg_regions,  ///< TODOCUMENT
                                                                              const ostream_ref_opt    &arg_ostream   ///< An optional reference to an ostream to which any logging should be sent
                                                                              ) {
	if ( arg_pdb_list.size() != arg_regions.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Number of regions lists must match number of PDBs in pdb_list to restrict."));
	}
	pdb_list new_pdb_list;
	new_pdb_list.reserve( arg_pdb_list.size() );

	// \TODO Come C++17 and structured bindings, use here
	for (const boost::tuple<const region_vec_opt &, const pdb &> &the_pair : combine( arg_regions, arg_pdb_list ) ) {
		new_pdb_list.push_back(
			backbone_complete_subset_of_pdb(
				get_regions_limited_pdb(
					the_pair.get<0>(),
					the_pair.get<1>()
				),
				arg_ostream
			).first
		);
	}
	return new_pdb_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
///
/// \relates protein_list
///
/// \TODO Consider taking an ostream_ref_opt argument rather than assuming cerr
///       (fix all errors, *then* provide default of boost::none)
protein_list cath::file::build_protein_list_of_pdb_list(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                        ) {
	protein_list new_protein_list;
	new_protein_list.reserve(arg_pdb_list.size());
	for (const pdb &the_pdb : arg_pdb_list) {
		new_protein_list.push_back( build_protein_of_pdb( the_pdb, ref( cerr ) ).first );
	}
	return new_protein_list;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
///
/// \relates protein_list
protein_list cath::file::build_protein_list_of_pdb_list_and_names(const pdb_list      &arg_pdb_list, ///< TODOCUMENT
                                                                  const name_set_list &arg_name_sets ///< TODOCUMENT
                                                                  ) {
	const size_t num_names = arg_name_sets.size();
	if ( arg_pdb_list.size() != num_names ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to make proteins from pdb_list and names because the numbers mismatch"));
	}
	protein_list new_proteins = build_protein_list_of_pdb_list( arg_pdb_list );
	for (boost::tuple<protein &, const name_set &> &&x : combine( new_proteins, arg_name_sets ) ) {
		x.get<0>().set_name_set( x.get<1>() );
	}
	return new_proteins;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
amino_acid_vec_vec cath::file::get_amino_acid_lists(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                    ) {
	amino_acid_vec_vec amino_acid_lists;
	amino_acid_lists.reserve( arg_pdb_list.size() );
	for (const pdb &the_pdb : arg_pdb_list) {
		amino_acid_lists.push_back( get_amino_acid_list( the_pdb ) );
	}
	return amino_acid_lists;
}

/// \brief TODOCUMENT
///
/// \relates pdb_list
residue_id_vec_vec cath::file::get_residue_ids_of_first_chains__backbone_unchecked(const pdb_list &arg_pdb_list ///< TODOCUMENT
                                                                                   ) {
	const size_t num_pdbs = arg_pdb_list.size();
	residue_id_vec_vec residue_ids;
	residue_ids.reserve( num_pdbs );
	for (const size_t &pdb_ctr : indices( num_pdbs ) ) {
		residue_ids.push_back( arg_pdb_list[ pdb_ctr ].get_residue_ids_of_first_chain__backbone_unchecked() );
	}
	return residue_ids;
}

/// \brief Get a cache of the indices of the backbone-complete indices for the specified PDBs
///
/// \relates pdb_list
backbone_complete_indices_vec cath::file::get_backbone_complete_indices(const pdb_list &arg_pdbs ///< The PDBs to query
                                                                        ) {
	return transform_build<backbone_complete_indices_vec>(
		arg_pdbs,
		[] (const pdb &x) { return get_backbone_complete_indices( x ); }
	);
}

/// \brief Get the lists of residue IDs for all chains of each of the PDBs in the specified pdb_list
///
/// \relates pdb_list
residue_id_vec_vec cath::file::get_backbone_complete_residue_ids(const pdb_list &arg_pdb_list ///< The PDBs to query
                                                                 ) {
	return transform_build<residue_id_vec_vec>(
		indices( arg_pdb_list.size() ),
		[&] (const size_t &x) {
			return get_backbone_complete_residue_ids( arg_pdb_list[ x ] );
		}
	);
}

/// \brief Get the lists of residue IDs for the first chains of each of the PDBs in the specified pdb_list
///
/// \relates pdb_list
residue_id_vec_vec cath::file::get_backbone_complete_residue_ids_of_first_chains(const pdb_list &arg_pdb_list ///< The PDBs to query
                                                                                 ) {
	return transform_build<residue_id_vec_vec>(
		indices( arg_pdb_list.size() ),
		[&] (const size_t &x) {
			return get_backbone_complete_residue_ids_of_first_chain( arg_pdb_list[ x ] );
		}
	);
}
