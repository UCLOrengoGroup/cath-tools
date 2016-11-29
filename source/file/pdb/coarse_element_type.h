/// \file
/// \brief The coarse_element_type class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_PDB_COARSE_ELEMENT_TYPE_H
#define _CATH_TOOLS_SOURCE_FILE_PDB_COARSE_ELEMENT_TYPE_H

#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief Coarsely represent an element type as one of the PDB core element types (C/CA/CB/N/O) or NON_CORE
		enum class coarse_element_type : char {
			CARBON,       ///< Represent a carbon atom
			CARBON_ALPHA, ///< Represent a carbon_alpha atom
			CARBON_BETA,  ///< Represent a carbon_beta atom
			NITROGEN,     ///< Represent a nitrogen atom
			OXYGEN,       ///< Represent a oxygen atom
			NON_CORE      ///< Represent a non-core atom
		};

		std::string to_string(const coarse_element_type &);

		std::ostream & operator<<(std::ostream &,
		                          const coarse_element_type &);

	} // namespace file
} // namespace cath

#endif
