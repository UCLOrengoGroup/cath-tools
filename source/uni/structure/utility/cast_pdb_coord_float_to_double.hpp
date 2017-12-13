/// \file
/// \brief The cast_pdb_coord_float_to_double header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_UTILITY_CAST_PDB_COORD_FLOAT_TO_DOUBLE_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_UTILITY_CAST_PDB_COORD_FLOAT_TO_DOUBLE_HPP

namespace cath {

	/// \brief Convert a float representing PDB coord dimension value into a double
	///
	/// Motivation: a simple static_cast won't always be enough to remove the
	/// rounding error introduced when the value was converted to a float.
	///
	/// eg :
	/// `static_cast<double>           ( 0.1f )` is 0.100000001490116119384765625
	/// `cast_pdb_coord_float_to_double( 0.1f )` is 0.100000000000000005551115123
	inline double cast_pdb_coord_float_to_double(const float &arg) {
		return static_cast<double>( static_cast<int>( 1000.0f * arg + 0.5 ) ) / 1000.0;
	}

} // namespace cath

#endif
