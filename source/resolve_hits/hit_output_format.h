/// \file
/// \brief The hit class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_HIT_OUTPUT_FORMAT_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_HIT_OUTPUT_FORMAT_H

namespace cath {
	namespace rslv {

		/// \brief Represent the different formats in which hits can be output
		enum class hit_output_format {
			CLASS, ///< Default format for outputting C++ classes
			JON    ///< Format in which hits should be output for Jon's use
		};


	} // namespace rslv
} // namespace cath

#endif
