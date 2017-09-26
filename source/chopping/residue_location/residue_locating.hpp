/// \file
/// \brief The residue_locating header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATING_H
#define _CATH_TOOLS_SOURCE_CHOPPING_RESIDUE_LOCATION_RESIDUE_LOCATING_H

#include <iosfwd>
#include <string>

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		enum class residue_locating : char {
			INDEX,         ///< TODOCUMENT
			NAME,          ///< TODOCUMENT
			NAME_AND_INDEX ///< TODOCUMENT
		};

		std::string to_string(const residue_locating &);
		std::ostream & operator<<(std::ostream &,
		                          const residue_locating &);

		residue_locating make_residue_locating_of_has_name_and_has_index(const bool &,
		                                                                 const bool &);

	} // namespace chop
} // namespace cath

#endif
