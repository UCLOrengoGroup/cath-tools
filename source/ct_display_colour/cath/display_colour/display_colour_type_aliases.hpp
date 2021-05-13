/// \file
/// \brief The display_colour type aliases header

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

#ifndef _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_TYPE_ALIASES_HPP

#include <map>
#include <optional>
#include <vector>

#include <boost/config.hpp> /// \todo Come a resolution for Boost Trac tickets 12142 & 12179, remove this #include

#include "cath/common/type_aliases.hpp"

namespace cath {
	class display_colour;

	/// \brief TODOCUMENT
	using display_colour_opt = ::std::optional<display_colour>;

	/// \brief TODOCUMENT
	using display_colour_set = std::set<display_colour>;

	/// \brief TODOCUMENT
	using display_colour_vec = std::vector<display_colour>;


    /// \brief TODOCUMENT
	using size_display_colour_map = std::map<size_t, display_colour>;
	/// \brief TODOCUMENT
	using size_display_colour_map_val = size_display_colour_map::value_type;


    /// \brief TODOCUMENT
	using size_size_display_colour_map = std::map<size_size_pair, display_colour>;
	/// \brief TODOCUMENT
	using size_size_display_colour_map_val = size_size_display_colour_map::value_type;

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_TYPE_ALIASES_HPP
