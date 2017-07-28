/// \file
/// \brief The restrict_to_single_linkage_extension class definitions

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

#include "restrict_to_single_linkage_extension.hpp"

#include "common/cpp14/cbegin_cend.hpp"
#include "structure/geometry/coord.hpp"

using namespace cath;
using namespace cath::geom;

/// \brief Find the coords within the specified vector that can be connected by a chain of links of
///        the specified distance or shorter (via the other coords) to one of the coords before the
///        specified point. Move all those coords (and the original coords) to the start of the vector
///        and erase the rest.
///
/// More informally:
/// * I'll give you a vector of coords, with a core of coords at the start indicated by an
///   iterator at the end and a distance.
/// * You keep extending the core to include any coords within the distance
///   of the current core. Modify the vector to contain only those coords.
///
/// The coords may be in any order in the final vector.
void cath::geom::restrict_to_single_linkage_extension(coord_vec           &arg_coords,            ///< The vector to modify
                                                      const coord_vec_itr &arg_core_end_itr,      ///< The point in arg_coords at the end of the core coords and at the start of the non-core coords
                                                      const double        &arg_extension_distance ///< The maximum distance by which the core may be extended in each step
                                                      ) {
	const double distance_sq     = arg_extension_distance * arg_extension_distance;
	const auto   end_itr         = std::end  ( arg_coords );
	auto         added_begin_itr = std::begin( arg_coords );
	auto         added_end_itr   = arg_core_end_itr;

	// Keep looping while some were added in the previous iteration (or, on the first pass, if there is a non-empty core)
	// and there are more to add
	while ( added_begin_itr != added_end_itr && added_end_itr != common::cend( arg_coords ) ) {
		// Partition the coords outside the core according to whether they're within distance
		// of any of the most-recently added coords
		// (ie move all those coords within distance to the start of the range)
		const auto new_added_end_itr = std::partition(
			added_end_itr,
			end_itr,
			[&] (const coord &candidate) {
				return std::any_of(
					added_begin_itr,
					added_end_itr,
					[&] (const coord &recently_added) {
						return squared_distance_between_points( candidate, recently_added ) <= distance_sq;
					}
				);
			}
		);
		// Updated the iterators at the begin/end of the most-recently added range
		added_begin_itr = added_end_itr;
		added_end_itr   = new_added_end_itr;
	}

	// Remove any coords that never got added
	arg_coords.erase(
		added_end_itr,
		common::cend( arg_coords )
	);
}

/// \brief Find the coords within the specified vector that can be connected by a chain of links of
///        the specified distance or shorter (via the other coords) to one of the coords before the
///        specified point. Return a vector that contains all those core coords and no others.
///
/// More informally:
/// * I'll give you a vector of coords, with a core of coords at the start indicated by an
///   iterator at the end and a distance.
/// * You keep extending the core to include any coords within the distance
///   of the current core. Modify the vector to contain only those coords.
///
/// The coords may be in any order in the final vector.
coord_vec cath::geom::restrict_to_single_linkage_extension_copy(coord_vec     arg_coords,            ///< The vector to copy and modify
                                                                const size_t &arg_core_end_offset,   ///< The index of the point in arg_coords at the end of the core coords and at the start of the non-core coords
                                                                const double &arg_extension_distance ///< The maximum distance by which the core may be extended in each step
                                                                ) {
	restrict_to_single_linkage_extension(
		arg_coords,
		std::next( begin( arg_coords ), static_cast<ptrdiff_t>( arg_core_end_offset ) ),
		arg_extension_distance
	);
	return arg_coords;
}
