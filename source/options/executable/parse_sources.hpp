/// \file
/// \brief The parse_sources class header

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

#ifndef _CATH_TOOLS_SOURCE_OPTIONS_EXECUTABLE_PARSE_SOURCES_H
#define _CATH_TOOLS_SOURCE_OPTIONS_EXECUTABLE_PARSE_SOURCES_H

namespace cath {
	namespace opts {

		/// \brief Represent the sources from which options should be parsed
		enum class parse_sources : bool {
			CMND_LINE_ONLY,   ///< Only parse options from the command-line
			CMND_ENV_AND_FILE ///< Parse options from the command-line, environment-variables and configuration file
		};

	} // namespace opts
} // namespace cath

#endif
