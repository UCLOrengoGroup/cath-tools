/// \file
/// \brief The data_option class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_OPTION_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_OPTION_HPP

#include <iosfwd>

namespace cath::opts::detail {

	/// \brief Represent the different options that a data_dirs_options_block supports for each file type
	enum class data_option : char {
		PATH,   ///< Enum value to denote a file type's path option
		PREFIX, ///< Enum value to denote a file type's prefix option
		SUFFIX  ///< Enum value to denote a file type's suffix option
	};

	std::ostream & operator<<(std::ostream &,
	                          const data_option &);

} // namespace cath::opts::detail

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_OPTIONS_DATA_OPTION_HPP
