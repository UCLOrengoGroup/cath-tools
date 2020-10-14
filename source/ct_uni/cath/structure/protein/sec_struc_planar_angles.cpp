/// \file
/// \brief The sec_struc_planar_angles class definitions

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

#include "sec_struc_planar_angles.hpp"

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/throw_exception.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/file/sec/sec_file_record.hpp"
#include "cath/structure/geometry/coord.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;
using ::std::string;

const sec_struc_planar_angles sec_struc_planar_angles::NULL_SEC_STRUC_PLANAR_ANGLES(0.0, 0.0, 0.0);

/// \brief Ctor for sec_struc_planar_angles
sec_struc_planar_angles::sec_struc_planar_angles(const double &prm_angle_x,       ///< The angle on the plane defined by the x-axis
                                                 const double &prm_angle_minus_y, ///< The angle on the plane defined by the negative y-axis
                                                 const double &prm_angle_z        ///< The angle on the plane defined by the z-axis
                                                 ) : planar_angle_x      (prm_angle_x       ),
                                                     planar_angle_minus_y(prm_angle_minus_y ),
                                                     planar_angle_z      (prm_angle_z       ) {
	using ::boost::math::isfinite;
	if (!isfinite(planar_angle_x) || !isfinite(planar_angle_minus_y) || !isfinite(planar_angle_z)) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Arguments angle_x, angle_y and angle_z must be a normal, finite floating-point numbers"));
	}
}

/// \brief Getter for the angle on the plane defined by the x-axis
double sec_struc_planar_angles::get_planar_angle_x() const {
	return planar_angle_x;
}

/// \brief Getter for the angle on the plane defined by the negative y-axis
double sec_struc_planar_angles::get_planar_angle_minus_y() const {
	return planar_angle_minus_y;
}

/// \brief Getter for the angle on the plane defined by the z-axis
double sec_struc_planar_angles::get_planar_angle_z() const {
	return planar_angle_z;
}

/// \brief Generate a string describing the specified sec_struc_planar_angles
///
/// \relates sec_struc_planar_angles
string cath::to_string(const sec_struc_planar_angles &prm_sec_struc_planar_angles ///< The sec_struc_planar_angles to describe
                       ) {
	return "sec_struc_planar_angles[ "
		+ std::to_string( prm_sec_struc_planar_angles.get_planar_angle_x()       )
		+ ", "
		+ std::to_string( prm_sec_struc_planar_angles.get_planar_angle_minus_y() )
		+ ", "
		+ std::to_string( prm_sec_struc_planar_angles.get_planar_angle_z()       )
		+ " ]";
}

/// \brief Calculate the planar angles between the two specified sec_file_records
///
/// \relates sec_struc_planar_angles
///
/// \relatesalso sec_file_record
sec_struc_planar_angles cath::make_planar_angles(const sec_file_record &prm_sec_1, ///< The first  sec_file_record
                                                 const sec_file_record &prm_sec_2  ///< The second sec_file_record
                                                 ) {
	return {
		angle_in_degrees( planar_angle_between( coord{  1,  0,  0 }, prm_sec_1.get_unit_dirn(), prm_sec_2.get_unit_dirn() ) ),
		angle_in_degrees( planar_angle_between( coord{  0, -1,  0 }, prm_sec_1.get_unit_dirn(), prm_sec_2.get_unit_dirn() ) ),
		angle_in_degrees( planar_angle_between( coord{  0,  0,  1 }, prm_sec_1.get_unit_dirn(), prm_sec_2.get_unit_dirn() ) )
	};
}
