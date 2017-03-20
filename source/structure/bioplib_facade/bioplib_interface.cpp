/// \file
/// \brief The bioplib interface definitions

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

#include "bioplib_interface.hpp"

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/irange.hpp>

#include "common/size_t_literal.hpp"
#include "exception/runtime_error_exception.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/geometry/rotation.hpp"

extern "C" {
	#include "bioplib/fit.h"
}

using namespace cath;
using namespace cath::common;
using namespace cath::geom;
using namespace std;

using boost::irange;
using boost::numeric_cast;

/// \brief Access to bioplib's matfit() superposing subroutine as used by cath::superposition
///
/// \relates bioplib_coord_list_interfacer
///
/// Rather than using this subroutine directly, access it via cath::superposition instead
rotation cath::bioplib_fit(const coord_list &arg_coord_list_1, ///< TODOCUMENT
                           const coord_list &arg_coord_list_2  ///< TODOCUMENT
                           ) {
	const size_t num_coords = check_non_empty_and_equal_size(
		arg_coord_list_1,
		arg_coord_list_2
	);

	// Build appropriate bioplib_coord_list_interfacers for each of the coord_lists
	bioplib_coord_list_interfacer coord_list_interfacer_1(arg_coord_list_1);
	bioplib_coord_list_interfacer coord_list_interfacer_2(arg_coord_list_2);

	// Calculate rotation matrix
	double temp_superposition_matrix[3][3];
	const bool matfit_return_value = numeric_cast<bool>(matfit(
		coord_list_interfacer_1.get_interface(),
		coord_list_interfacer_2.get_interface(),
		temp_superposition_matrix,
		numeric_cast<int>(num_coords),
		nullptr,
		numeric_cast<short int>(false)
	));

	doub_vec untransposed_rotation_matrix;
	untransposed_rotation_matrix.reserve(coord::NUM_DIMS * coord::NUM_DIMS);
	for (const size_t &ctr1 : irange( 0_z, coord::NUM_DIMS ) ) {
		for (const size_t &ctr2 : irange( 0_z, coord::NUM_DIMS ) ) {
			// NOTE: Switching indices here because bioplib appears to index in the opposite way
			untransposed_rotation_matrix.push_back(temp_superposition_matrix[ctr2][ctr1]);
		}
	}
	const rotation sup_rotation(untransposed_rotation_matrix);

	// If Rotation matrix doesn't work for some reason, return -1
	// \todo Test this error by sending a 1 point list
	if (!matfit_return_value) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Failure to superpose using call to ProFit's/bioplib's fit()"));
	}

	return sup_rotation;
}

/// \brief TODOCUMENT
bioplib_coord_list_interfacer::bioplib_coord_list_interfacer(const coord_list &arg_coord_list ///< TODOCUMENT
                                                             ) {
	coors.reserve(arg_coord_list.size());
	for (const coord &loop_coord : arg_coord_list) {
		coors.push_back(COOR_of_coord(loop_coord));
	}
}

/// \brief TODOCUMENT
COOR * bioplib_coord_list_interfacer::get_interface() {
	return &coors.front();
}

/// \brief TODOCUMENT
///
/// \relates bioplib_coord_list_interfacer
COOR cath::COOR_of_coord(const coord &arg_coord ///< TODOCUMENT
                         ) {
	COOR result;
	result.x = arg_coord.get_x();
	result.y = arg_coord.get_y();
	result.z = arg_coord.get_z();
	return result;
}

/// \brief TODOCUMENT
///
/// \relates bioplib_coord_list_interfacer
coord cath::coord_of_COOR(const COOR &arg_coor ///< TODOCUMENT
                          ) {
	return coord(arg_coor.x, arg_coor.y, arg_coor.z);
}

