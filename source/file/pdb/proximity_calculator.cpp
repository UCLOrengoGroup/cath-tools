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

#include "common/algorithm/copy_build.hpp"
#include "chopping/region/regions_limiter.hpp"
#include "common/algorithm/sort_uniq_build.hpp"
#include "common/boost_addenda/range/front.hpp"
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
using boost::irange;
using std::tie;

/// \brief Get the res_index_key_coord_and_dist values for the specified PDB, restricted
///        to the specified regions
res_index_key_coord_and_dist_vec proximity_calculator::get_res_indices_key_coords_and_dists(const pdb            &arg_pdb,    ///< The PDB to index
                                                                                            const region_vec_opt &arg_regions ///< The regions of the PDB to index
                                                                                            ) {
	detail::res_index_key_coord_and_dist_vec results;
	regions_limiter the_limiter{ arg_regions };

	// Loop over the residues (tracking the index)
	for (const size_t &res_ctr : irange( 0_z, arg_pdb.get_num_residues() ) ) {
		const auto &res        = arg_pdb.get_residue_of_index__backbone_unchecked( res_ctr );

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
proximity_calculator::proximity_calculator(const pdb            &arg_pdb,    ///< The PDB to index
                                           const region_vec_opt &arg_regions ///< The regions to index
                                           ) : source_pdb { arg_pdb },
                                               obj_res_indices_key_coords_and_dists{
                                               	get_res_indices_key_coords_and_dists( arg_pdb, arg_regions )
                                               } {
}

/// \brief Return whether there are any atoms in the original PDB that are
///        within the specified distance of the specified coord
bool proximity_calculator::is_within_distance(const coord  &arg_coord,   ///< The coord to query
                                              const double &arg_distance ///< The maximum distance to the original PDB
                                              ) const {
	// Prepare the squared distance for quicker calculation later
	const double distance_sq = arg_distance * arg_distance;

	// Return whether the lambda returns true for any of the obj_res_indices_key_coords_and_dists
	return any_of(
		obj_res_indices_key_coords_and_dists,
		[&] (const res_index_key_coord_and_dist &res) {

			// Prepare some information for the residue
			const size_t &res_indx     = res.res_index;
			const coord  &res_locn     = res.key_coord;
			const double &res_dist     = res.furthest_atom_dist;
			const auto    squared_dist = squared_distance_between_points( res_locn, arg_coord );

			// If the distance to the key coord is less than the distance then the condition
			// is met, so return true
			if ( squared_dist <= distance_sq ) {
				return true;
			}

			// If the distance to the key coord is longer than the allowable distance plus the
			// furthest distance from any of the atoms to the key coord, then none of the atoms
			// can possibly be close enough, so return false
			if ( squared_dist > ( arg_distance + res_dist ) * ( arg_distance + res_dist ) ) {
				return false;
			}

			// Otherwise, return whether any of the atoms are close enough
			return any_of(
				source_pdb.get().get_residue_of_index__backbone_unchecked( res_indx ),
				[&] (const pdb_atom &atom) {
					return ( squared_distance_between_points( atom.get_coord(), arg_coord ) <= distance_sq );
				}
			);
		}
	);
}

/// \brief Return whether there are any atoms in the original PDB that are
///        within the specified distance of the specified pdb_residue
///
/// \relates proximity_calculator
bool cath::file::is_within_distance(const proximity_calculator &arg_prox_calc, ///< The proximity_calculator to query
                                    const pdb_residue          &arg_residue,   ///< The pdb_residue to query
                                    const double               &arg_distance   ///< The maximum distance to the original PDB
                                    ) {
	return any_of(
		arg_residue,
		[&] (const pdb_atom &x) {
			return arg_prox_calc.is_within_distance( x.get_coord(), arg_distance );
		}
	);
}


/// \brief Return the labels of any DNA/RNA chains in the specified PDB that come within
///        the specified distance of the original PDB of the specified proximity_calculator
///
/// \relates proximity_calculator
chain_label_set cath::file::nearby_dna_rna_chain_labels(const pdb                  &arg_pdb,       ///< The PDB in question
                                                        const proximity_calculator &arg_prox_calc, ///< The proximity_calculator to query
                                                        const double               &arg_distance   ///< The maximum distance to the original PDB
                                                        ) {
	chain_label_set results;

	// Loop over residues
	for (const pdb_residue &res : arg_pdb) {

		// If the residue is a DNA/RNA residue on a not-yet-recorded chain
		if ( res.get_amino_acid().get_type() == amino_acid_type::DNA && ! contains( results, get_chain_label( res ) ) ) {

			// If the residue has any atoms within distance of the arg_prox_calc
			if ( is_within_distance( arg_prox_calc, res, arg_distance ) ) {

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
chain_label_set cath::file::nearby_dna_rna_chain_labels(const pdb                  &arg_pdb,       ///< The PDB in question
                                                        const proximity_calculator &arg_prox_calc, ///< The proximity_calculator to query
                                                        const doub_opt             &arg_dist_opt   ///< The maximum distance to the original PDB
                                                        ) {
	if ( arg_dist_opt ) {
		return nearby_dna_rna_chain_labels( arg_pdb, arg_prox_calc, arg_dist_opt );
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
void cath::file::restrict_to_linkage_proximate(coord_vec                  &arg_coords,            ///< The coord vector to modify
                                               const proximity_calculator &arg_prox_calc,         ///< The proximity_calculator to query
                                               const double               &arg_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                               const double               &arg_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                               ) {
	restrict_to_single_linkage_extension(
		arg_coords,
		[&] () {
			return partition(
				arg_coords,
				[&] (const coord &x) {
					return arg_prox_calc.is_within_distance(
						x,
						arg_primary_distance
					);
				}
			);
		} (),
		arg_extension_distance
	);
}

/// \brief Return a copy of the specified coord_vec, restricted to those coords that are within the specified primary distance
///        to the original PDB of the specified proximity_calculator and those coords that can be reached
///        from them through other coords in steps of the specified extension distance or less
///
/// The returned data may be unsorted
///
/// \relates proximity_calculator
coord_vec cath::file::restrict_to_linkage_proximate_copy(coord_vec                   arg_coords,            ///< The coord vector to copy and return modified version of
                                                         const proximity_calculator &arg_prox_calc,         ///< The proximity_calculator to query
                                                         const double               &arg_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                                         const double               &arg_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                                         ) {
	restrict_to_linkage_proximate( arg_coords, arg_prox_calc, arg_primary_distance, arg_extension_distance );
	return arg_coords;
}

/// \brief Get a pdb_residue_vec of those post-TER residues in the specified PDB that are
///        nearby to the original PDB of the specified proximity_calculator (within the specified primary distance
///        or reachable from there through other coords in steps of the specified extension distance or less)
pdb_residue_vec cath::file::get_nearby_post_ter_res_atoms(const pdb                  &arg_pdb,               ///< The PDB with the post-TER records to query
                                                          const proximity_calculator &arg_prox_calc,         ///< The proximity_calculator to query
                                                          const double               &arg_primary_distance,  ///< The maximum primary allowed distance to the PDB
                                                          const double               &arg_extension_distance ///< The maximum extension distance to pull in more coords in a single-linkage fashion
                                                          ) {
	// Create a function object for less-than-comparing coords
	//
	// (this isn't an operator<() for coord because it doesn't make sense in most
	//  contexts but it makes sense for sorting a list of coords so that the list can
	//  be binary_search()-ed)
	const auto coord_less_fn = [] (const coord &lhs, const coord &rhs) {
		return (
			tie( lhs.get_x(), lhs.get_y(), lhs.get_z() )
			<
			tie( rhs.get_x(), rhs.get_y(), rhs.get_z() )
		);
	};

	// Get a sorted list of the post-ter residue coords that are "nearby"
	// (within arg_primary_distance) according to proximity_calculator or
	// are connected to any such through a chain of others with a maximum
	// link distance of arg_extension_distance
	const coord_vec resorted_nearbys = sort_copy(
		restrict_to_linkage_proximate_copy(
			get_all_coords( arg_pdb.get_post_ter_residues() ),
			arg_prox_calc,
			arg_primary_distance,
			arg_extension_distance
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
						atom.get_coord(),
						coord_less_fn
					);
				} )
			)
		};
	};

	const auto res_is_not_empty_fn = [] (const pdb_residue &res) { return ! res.empty(); };

	// Return the ???
	return copy_build<pdb_residue_vec>(
		arg_pdb.get_post_ter_residues()
			| transformed( restrict_res_to_nearby_atoms_fn )
			| filtered   ( res_is_not_empty_fn             )
	);
}


