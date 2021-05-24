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

// clang-format off
namespace cath::file { class sec_file_record; }
// clang-format on

#include <string>

#include "cath/common/exception/invalid_argument_exception.hpp"

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
		constexpr sec_struc_planar_angles( const double &, const double &, const double & );
		[[nodiscard]] constexpr double get_planar_angle_x() const;
		[[nodiscard]] constexpr double get_planar_angle_minus_y() const;
		[[nodiscard]] constexpr double get_planar_angle_z() const;
	};

	/// \brief Ctor for sec_struc_planar_angles
	///
	/// \param prm_angle_x       The angle on the plane defined by the x-axis
	/// \param prm_angle_minus_y The angle on the plane defined by the negative y-axis
	/// \param prm_angle_z       The angle on the plane defined by the z-axis
	constexpr sec_struc_planar_angles::sec_struc_planar_angles( const double &prm_angle_x,
	                                                            const double &prm_angle_minus_y,
	                                                            const double &prm_angle_z ) :
	        planar_angle_x( prm_angle_x ), planar_angle_minus_y( prm_angle_minus_y ), planar_angle_z( prm_angle_z ) {
		if ( !( planar_angle_x > ::std::numeric_limits<double>::lowest()
		        && planar_angle_x < ::std::numeric_limits<double>::max() )
		     || !( planar_angle_minus_y > ::std::numeric_limits<double>::lowest()
		           && planar_angle_minus_y < ::std::numeric_limits<double>::max() )
		     || !( planar_angle_z > ::std::numeric_limits<double>::lowest()
		           && planar_angle_z < ::std::numeric_limits<double>::max() ) ) {
			BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
			  "Arguments angle_x, angle_y and angle_z must be a normal, finite floating-point numbers" ) );
		}
	}

	/// \brief Getter for the angle on the plane defined by the x-axis
	constexpr double sec_struc_planar_angles::get_planar_angle_x() const {
		return planar_angle_x;
	}

	/// \brief Getter for the angle on the plane defined by the negative y-axis
	constexpr double sec_struc_planar_angles::get_planar_angle_minus_y() const {
		return planar_angle_minus_y;
	}

	/// \brief Getter for the angle on the plane defined by the z-axis
	constexpr double sec_struc_planar_angles::get_planar_angle_z() const {
		return planar_angle_z;
	}

	inline constexpr sec_struc_planar_angles NULL_SEC_STRUC_PLANAR_ANGLES( 0.0, 0.0, 0.0 );

	std::string to_string(const sec_struc_planar_angles &);

	sec_struc_planar_angles make_planar_angles(const file::sec_file_record &,
	                                           const file::sec_file_record &);

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_SEC_STRUC_PLANAR_ANGLES_HPP
