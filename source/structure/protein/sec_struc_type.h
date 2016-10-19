/// \file
/// \brief The sec_struc_type header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_SEC_STRUC_TYPE_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_PROTEIN_SEC_STRUC_TYPE_H

#include <iosfwd>

namespace cath {

	/// \brief Represent the possible secondary structure types
	enum class sec_struc_type {
		ALPHA_HELIX,
		BETA_STRAND,
		COIL
	};

	std::istream & operator>>(std::istream &,
	                          sec_struc_type &);

	std::ostream & operator<<(std::ostream &,
	                          const sec_struc_type &);

} // namespace cath

#endif
