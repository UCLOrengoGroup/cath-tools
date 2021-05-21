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

#include <boost/throw_exception.hpp>

#include "cath/file/sec/sec_file_record.hpp"
#include "cath/structure/geometry/coord.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::geom;
using ::std::string;

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
