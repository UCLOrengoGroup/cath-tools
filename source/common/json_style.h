/// \file
/// \brief The json_style class header

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

#ifndef _CATH_TOOLS_SOURCE_OUTPUTTER_SUPERPOSITION_OUTPUTTER_SUP_JSON_STYLE_H
#define _CATH_TOOLS_SOURCE_OUTPUTTER_SUPERPOSITION_OUTPUTTER_SUP_JSON_STYLE_H

namespace cath {
	namespace common {

		/// \brief The style in which superposition JSON should be written
		///
		/// \todo This should be used in more superposition/JSON code rather than just in json_style.
		///       Extend the use of this all the way to the call to the Boost code.
		enum class json_style : bool {
			PRETTY, ///< Insert white-space characters to format the JSON in a more human-readable layout
			COMPACT ///< Suppress superfluous white-space characters to keep the JSON compact
		};

	} // namespace common
} // namespace cath

#endif
