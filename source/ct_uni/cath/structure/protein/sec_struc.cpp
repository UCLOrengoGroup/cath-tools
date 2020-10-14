/// \file
/// \brief The sec_struc class definitions

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

#include "sec_struc.hpp"

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/operators.hpp>
#include <boost/throw_exception.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::geom;
using namespace ::std;

///// \brief A static black sec_struc object
//const sec_struc sec_struc::NULL_SEC_STRUC(
//	0,
//	0,
//	sec_struc_type::ALPHA_HELIX,
//	ORIGIN_COORD,
//	ORIGIN_COORD
//);

/// \brief TODOCUMENT
void sec_struc::check_planar_angles_index_is_valid(const size_t &prm_index ///< TODOCUMENT
                                                   ) const {
	if (prm_index >= get_num_planar_angles()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Pair details index is out of range"));
	}
}

/// \brief TODOCUMENT
void sec_struc::check_sec_struc_type() const {
	if (type == sec_struc_type::COIL) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("sec_struc's type cannot be sec_struc_type::COIL"));
	}
}

/// \brief TODOCUMENT
sec_struc::sec_struc(const size_t         &prm_start_residue_num, ///< TODOCUMENT
                     const size_t         &prm_stop_residue_num,  ///< TODOCUMENT
                     const sec_struc_type &prm_sec_struc_type,    ///< TODOCUMENT
                     coord                 prm_midpoint,          ///< TODOCUMENT
                     coord                 prm_unit_dirn          ///< TODOCUMENT
                     ) : start_residue_num ( prm_start_residue_num      ),
                         stop_residue_num  ( prm_stop_residue_num       ),
                         type              ( prm_sec_struc_type         ),
                         midpoint          ( std::move( prm_midpoint  ) ),
                         unit_dirn         ( std::move( prm_unit_dirn ) ) {
	check_sec_struc_type();
}

/// \brief TODOCUMENT
void sec_struc::set_planar_angles(const sec_struc_planar_angles_vec &prm_planar_angles ///< TODOCUMENT
                                  ) {
	planar_angles.reserve(prm_planar_angles.size());
	planar_angles = prm_planar_angles;
}

/// \brief TODOCUMENT
size_t sec_struc::get_start_residue_num() const {
	return start_residue_num;
}

/// \brief TODOCUMENT
size_t sec_struc::get_stop_residue_num() const {
	return stop_residue_num;
}

/// \brief TODOCUMENT
sec_struc_type sec_struc::get_type() const {
	return type;
}

/// \brief TODOCUMENT
coord sec_struc::get_midpoint() const {
	return midpoint;
}

/// \brief TODOCUMENT
coord sec_struc::get_unit_dirn() const {
	return unit_dirn;
}

/// \brief TODOCUMENT
const sec_struc_planar_angles & sec_struc::get_planar_angles_of_index(const size_t &prm_index ///< TODOCUMENT
                                                                      ) const {
	check_planar_angles_index_is_valid(prm_index);
	return planar_angles[prm_index];
}

size_t sec_struc::get_num_planar_angles() const {
	return planar_angles.size();
}

/// \brief TODOCUMENT
///
/// \relates sec_struc
coord cath::calculate_inter_sec_struc_vector(const sec_struc &prm_src_sec_struc,   ///< The source secondary structure
                                             const sec_struc &prm_dest_sec_struc,  ///< The destination secondary structure
                                             const sec_struc &prm_anchor_sec_struc ///< The anchor secondary structure
                                             ) {
	// Calculate scalar/vector distances between secondary structures
	const rotation &anchor_frame = rotation_to_x_axis_and_x_y_plane(
		prm_src_sec_struc.get_unit_dirn(),
		prm_anchor_sec_struc.get_midpoint() - prm_src_sec_struc.get_midpoint()
	);
	return rotate_copy(
		anchor_frame,
		prm_dest_sec_struc.get_midpoint()   - prm_src_sec_struc.get_midpoint()
	);
}


