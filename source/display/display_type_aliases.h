/// \file
/// \brief The  class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef DISPLAY_TYPE_ALIASES_H_INCLUDED
#define DISPLAY_TYPE_ALIASES_H_INCLUDED

#include <boost/optional/optional_fwd.hpp>

#include "common/type_aliases.h"

#include <map>
#include <vector>

namespace cath {
	class display_colour;

	/// \brief TODOCUMENT
	using opt_display_colour = boost::optional<display_colour>;

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

}

#endif
