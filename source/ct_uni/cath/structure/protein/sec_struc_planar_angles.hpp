/// \file
/// \brief The sec_struc_planar_angles class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_PLANAR_ANGLES_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_PLANAR_ANGLES_HPP

namespace cath { namespace file { class sec_file_record; } }

#include <string>

namespace cath {

	// Planar angles between secondary structures as used in sec files
	class sec_struc_planar_angles final {
	private:
		/// \brief The y-z planar angle from the vector through the source sec_struc to the axis through the destination sec_struc
		double planar_angle_x;

		/// \brief The x-z planar angle from the vector through the source sec_struc to the axis through the destination sec_struc
		///
		/// Note that this is on the x-z plane (defined by -y) rather than the z-x plane (defined by y) as you might expect
		/// because secmake generates it this way to replicate previous behaviour.
		///
		/// \todo Check that this makes no difference because this is only be used for differences anyway.
		double planar_angle_minus_y;

		/// \brief The x-y planar angle from the vector through the source sec_struc to the axis through the destination sec_struc
		double planar_angle_z;

	public:
		sec_struc_planar_angles(const double &,
		                        const double &,
		                        const double &);
		[[nodiscard]] double get_planar_angle_x() const;
		[[nodiscard]] double get_planar_angle_minus_y() const;
		[[nodiscard]] double get_planar_angle_z() const;

		static const sec_struc_planar_angles NULL_SEC_STRUC_PLANAR_ANGLES;
	};

	std::string to_string(const sec_struc_planar_angles &);

	sec_struc_planar_angles make_planar_angles(const file::sec_file_record &,
	                                           const file::sec_file_record &);

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_PLANAR_ANGLES_HPP
