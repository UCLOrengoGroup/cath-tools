/// \file
/// \brief The proximity_calculator class definitions

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

#include "proximity_calculator.hpp"

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/partition.hpp>

#include "chopping/region/regions_limiter.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_build.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "file/pdb/pdb.hpp"
#include "structure/geometry/restrict_to_single_linkage_extension.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file::detail;
using namespace cath::file;
using namespace cath::geom;

using boost::adaptors::filtered;
using boost::adaptors::transformed;
using boost::algorithm::any_of;
using boost::range::binary_search;
using boost::range::partition;
using std::tie;

/// \brief Get the res_index_key_coord_and_dist values for the specified PDB, restricted
///        to the specified regions
res_index_key_coord_and_dist_vec proximity_calculator::get_res_indices_key_coords_and_dists(const pdb            &prm_pdb,    ///< The PDB to index
                                                                                            const region_vec_opt &prm_regions ///< The regions of the PDB to index
                                                                                            ) {
	detail::res_index_key_coord_and_dist_vec results;
	regions_limiter the_limiter{ prm_regions };

	// Loop over the residues (tracking the index)
	for (const size_t &res_ctr : indices( prm_pdb.get_num_residues() ) ) {
		const auto &res        = prm_pdb.get_residue_of_index__backbone_unchecked( res_ctr );

		// If this residue is included in the regions and isn't empty...
		const bool is_included = the_limiter.update_residue_is_included( res.get_residue_id() );
		if ( is_included && ! res.empty() ) {

			// Then grab the coordinate of a key atom (preferably the CA atom if it's present
			// or the first atom otherwise)
			const coord key_coord = res.has_carbon_alpha() ? get_carbon_alpha_coord( res )
			                                               : front( res ).get_coord();

			// Record this residue's index and key atom's coordinate, along with the maximum
			// distance from that coordinate to any of the other coordinates
			results.emplace_back(
				res_ctr,
				key_coord,
				max_dist_from_coord( res, key_coord )
			);
		}
	}
	return results;
}

/// \brief Ctor from a PDB and (optionally) some regions
proximity_calculator::proximity_calculator(const pdb            &prm_pdb,    ///< The PDB to index
                                           const region_vec_opt &prm_regions ///< The regions to index
                                           ) : source_pdb { prm_pdb },
                                               obj_res_indices_key_coords_and_dists{
                                               	get_res_indices_key_coords_and_dists( prm_pdb, prm_regions )
                                               } {
}

/// \brief Return whether there are any atoms in the original PDB that are
///        within the specified distance of the specified coord
bool proximity_calculator::is_within_distance(const coord  &prm_coord,   ///< The coord to query
                                              const double &prm_distance ///< The maximum distance to the original PDB
                                              ) const {
	// Prepare the squared distance for quicker calculation later
	const double distance_sq = prm_distance * prm_distance;

	// Return whether the lambda returns true for any of the obj_res_indices_key_coords_and_dists
	return any_of(
		obj_res_indices_key_coords_and_dists,
		[&] (const res_index_key_coord_and_dist &res) {

			// Prepare some information for the residue
			const size_t &res_indx     = res.res_index;
			const coord  &res_locn     = res.key_coord;
			const double &res_dist     = res.furthest_atom_dist;
			const auto    squared_dist = squared_distance_between_points( res_locn, prm_coord );

			// If the distance to the key coord is less than the distance then the condition
			// is met, so return true
			if ( squared_dist <= distance_sq ) {
				return true;
			}

			// If the distance to the key coord is longer than the allowable distance plus the
			// furthest distance from any of the atoms to the key coord, then none of the atoms
			// can possibly be close enough, so return false
			if ( squared_dist > ( prm_distance + res_dist ) * ( prm_distance + res_dist ) ) {
				return false;
			}

			// Otherwise, return whether any of the atoms are close enough
			return any_of(
				source_pdb.get().get_residue_of_index__backbone_unchecked( res_indx ),
				[&] (const pdb_atom &atom) {
					return ( squared_distance_between_points( atom.get_coord(), prm_coord ) <= distance_sq );
				}
			);
		}
	);
}

/// \brief Return whether there are any atoms in the original PDB that are
///        within the specified distance of the specified pdb_residue
///
/// \relates proximity_calculator
bool cath::file::is_within_distance(const proximity_calculator &prm_prox_calc, ///< The proximity_calculator to query
                                    const pdb_residue          &prm_residue,   ///< The pdb_residue to query
                                    const double               &prm_distance   ///< The maximum distance to the original PDB
                                    ) {
	return any_of(
		prm_residue,
		[&] (const pdb_atom &x) {
			return prm_prox_calc.is_within_distance( x.get_coord(), prm_distance );
		}
	);
}


/// \brief Return the labels of any DNA/RNA chains in the specified PDB that come within
///        the specified distance of the original PDB of the specified proximity_calculator
///
/// \relates proximity_calculator
chain_label_set cath::file::nearby_dna_rna_chain_labels(const pdb                  &prm_pdb,       ///< The PDB in question
                                                        const proximity_calculator &prm_prox_calc, ///< The proximity_calculator to query
                                                        const double               &prm_distance   ///< The maximum distance to the original PDB
                                                        ) {
	chain_label_set results;

	// Loop over residues
	for (const pdb_residue &res : prm_pdb) {

		// If the residue is a DNA/RNA residue on a not-yet-recorded chain
		if ( res.get_amino_acid().get_type() == amino_acid_type::DNA && ! contains( results, get_chain_label( res ) ) ) {

			// If the residue has any atoms within distance of the prm_prox_calc
			if ( is_within_distance( prm_prox_calc, res, prm_distance ) ) {

				// Insert the residue's chain_label into results
				results.insert( get_chain_label( res ) );
			}
		}
	}

	// Return results
	return results;
}

/// \brief Return the labels of any DNA/RNA chains in the specified PDB that come within
///        the specified distance of the original PDB of the specified proximity_calculator
///        or an empty set if the distance is none
///
/// \relates proximity_calculator
chain_label_set cath::file::nearby_dna_rna_chain_labels(const pdb                  &prm_pdb,       ///< The PDB in question
                                                        const proximity_calculator &prm_prox_calc, ///< The proximity_calculator to query
                                                        const doub_opt             &prm_dist_opt   ///< The maximum distance to the original PDB
                                                        ) {
	if ( prm_dist_opt ) {
		return nearby_dna_rna_chain_labels( prm_pdb, prm_prox_calc, *prm_dist_opt );
	}
	return {};
}

/// \brief Restrict the specified coord_vec to those coords that are within the specified primary distance
///        to the original PDB of the specified proximity_calculator and those coords that can be reached
///        from them through other coords in steps of the specified extension distance or less
///
/// The returned data may be unsorted
///
/// \relates proximity_calculator
void cath::file::restrict_to_linkage_proximate(coord_coord_linkage_pair_vec &prm_coords,            ///< The coord vector to modify
                                               const proximity_calculator   &prm_prox_calc,         ///< The proximity_calculator to query
                                               const double                 &prm_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                               const double                 &prm_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                               ) {
	restrict_to_single_linkage_extension(
		prm_coords,
		partition(
			prm_coords,
			[&] (const coord_coord_linkage_pair &x) {
				return prm_prox_calc.is_within_distance(
					x.first,
					prm_primary_distance
				);
			}
		),
		prm_extension_distance
	);
}

/// \brief Return a copy of the specified coord_coord_linkage_pair_vec, restricted to those coords that are
///        within the specified primary distance to the original PDB of the specified proximity_calculator
///        and those coords that can be reached from them through other coords in steps of the specified
///        extension distance or less
///
/// The returned data may be unsorted
///
/// \relates proximity_calculator
coord_coord_linkage_pair_vec cath::file::restrict_to_linkage_proximate_copy(coord_coord_linkage_pair_vec  prm_coords,            ///< The coord vector to copy and return modified version of
                                                                            const proximity_calculator   &prm_prox_calc,         ///< The proximity_calculator to query
                                                                            const double                 &prm_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                                                            const double                 &prm_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                                                            ) {
	restrict_to_linkage_proximate( prm_coords, prm_prox_calc, prm_primary_distance, prm_extension_distance );
	return prm_coords;
}

/// \brief Get a pdb_residue_vec of those post-TER residues in the specified PDB that are
///        nearby to the original PDB of the specified proximity_calculator (within the specified primary distance
///        or reachable from there through other coords in steps of the specified extension distance or less)
pdb_residue_vec cath::file::get_nearby_post_ter_res_atoms(const pdb                  &prm_pdb,               ///< The PDB with the post-TER records to query
                                                          const proximity_calculator &prm_prox_calc,         ///< The proximity_calculator to query
                                                          const double               &prm_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                                          const double               &prm_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                                          ) {
	// Create a function object for less-than-comparing the coord part of the coord_coord_linkage_pairs
	//
	// (this isn't an operator<() for coord because it doesn't make sense in most
	//  contexts but it makes sense for sorting a list of coords so that the list can
	//  be binary_search()-ed)
	const auto coord_less_fn = [] (const coord_coord_linkage_pair &lhs, const coord_coord_linkage_pair &rhs) {
		const auto &l = lhs.first;
		const auto &r = rhs.first;
		return (
			tie( l.get_x(), l.get_y(), l.get_z() )
			<
			tie( r.get_x(), r.get_y(), r.get_z() )
		);
	};

	// Get a sorted list of the post-ter residue coords that are "nearby"
	// (within prm_primary_distance) according to proximity_calculator or
	// are connected to any such through a chain of others with a maximum
	// link distance of prm_extension_distance
	const coord_coord_linkage_pair_vec resorted_nearbys = sort_copy(
		restrict_to_linkage_proximate_copy(
			get_all_coords_with_linkage( prm_pdb.get_post_ter_residues() ),
			prm_prox_calc,
			prm_primary_distance,
			prm_extension_distance
		),
		coord_less_fn
	);

	const auto restrict_res_to_nearby_atoms_fn = [&] (const pdb_residue &res) {
		return pdb_residue{
			res.get_residue_id(),
			copy_build<pdb_atom_vec>(
				res | filtered( [&] (const pdb_atom &atom) {
					return binary_search(
						resorted_nearbys,
						std::make_pair( atom.get_coord(), coord_linkage::ADD_ONLY ),
						coord_less_fn
					);
				} )
			)
		};
	};

	const auto res_is_not_empty_fn = [] (const pdb_residue &res) { return ! res.empty(); };

	// Return the ???
	return copy_build<pdb_residue_vec>(
		prm_pdb.get_post_ter_residues()
			| transformed( restrict_res_to_nearby_atoms_fn )
			| filtered   ( res_is_not_empty_fn             )
	);
}

/// \brief Get a pdb_residue_vec of those post-TER residues in the specified PDB that are
///        nearby to the original PDB of the specified proximity_calculator (within the specified primary distance
///        or reachable from there through other coords in steps of the specified extension distance or less)
///        or an empty set if the distance is none
///
/// \relates proximity_calculator
pdb_residue_vec cath::file::get_nearby_post_ter_res_atoms(const pdb                  &prm_pdb,               ///< The PDB with the post-TER records to query
                                                          const proximity_calculator &prm_prox_calc,         ///< The proximity_calculator to query
                                                          const doub_opt             &prm_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                                          const double               &prm_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                                          ) {
	if ( prm_primary_distance ) {
		return get_nearby_post_ter_res_atoms(
			prm_pdb,
			prm_prox_calc,
			*prm_primary_distance,
			prm_extension_distance
		);
	}
	return {};
}
