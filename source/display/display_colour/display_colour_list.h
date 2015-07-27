/// \file
/// \brief The display_colour_list class header

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

#ifndef DISPLAY_COLOUR_LIST_H_INCLUDED
#define DISPLAY_COLOUR_LIST_H_INCLUDED

#include "common/type_aliases.h"
#include "display/display_type_aliases.h"

#include <string>
#include <vector>

namespace cath {
	class display_colour;

	/// \brief TODOCUMENT
	class display_colour_list final {
	private:
		/// \brief TODOCUMENT
		display_colour_vec colours;

	public:
		explicit display_colour_list(const display_colour_vec &);

		size_t size() const;
		const display_colour & colour_of_index(const size_t &) const;

		static str_vec     DEFAULT_COLOURS_STRING_PARTS;
		static std::string COLOURS_SEPARATOR;
		static std::string DEFAULT_COLOURS_STRING;
	};

	const display_colour & colour_of_mod_index(const display_colour_list &,
	                                           const size_t &);

//	std::string name_of_colour_of_mod_index(const display_colour_list &,
//	                                        const display_colour_list::size_type &);

	display_colour_list make_display_colour_list_from_string(const std::string &);
}

#endif
