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

#include "sec_struc_planar_angles.h"

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/throw_exception.hpp>

#include "exception/invalid_argument_exception.h"

using namespace cath;
using namespace cath::common;

const sec_struc_planar_angles sec_struc_planar_angles::NULL_SEC_STRUC_PLANAR_ANGLES(0.0, 0.0, 0.0);

/// \brief Ctor for sec_struc_planar_angles.
sec_struc_planar_angles::sec_struc_planar_angles(const double &arg_angle_x,       ///< The angle on the plane defined by the x-axis
                                                 const double &arg_angle_minus_y, ///< The angle on the plane defined by the negative y-axis
                                                 const double &arg_angle_z        ///< The angle on the plane defined by the z-axis
                                                 ) : planar_angle_x      (arg_angle_x       ),
                                                     planar_angle_minus_y(arg_angle_minus_y ),
                                                     planar_angle_z      (arg_angle_z       ) {
	using boost::math::isfinite;
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
