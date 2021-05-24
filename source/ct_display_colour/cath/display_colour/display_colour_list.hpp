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

#ifndef _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_LIST_HPP
#define _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_LIST_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"

#include <string>
#include <vector>

// clang-format off
namespace cath { class display_colour; }
// clang-format on

namespace cath {

	/// \brief TODOCUMENT
	class display_colour_list final {
	private:
		/// \brief TODOCUMENT
		display_colour_vec colours;

	public:
		using const_iterator = display_colour_vec::const_iterator;

		explicit display_colour_list(display_colour_vec);

		[[nodiscard]] size_t                size() const;
		[[nodiscard]] const display_colour &colour_of_index( const size_t & ) const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	const display_colour & colour_of_mod_index(const display_colour_list &,
	                                           const size_t &);

	display_colour_list make_display_colour_list_from_string(const std::string &);

	display_colour_list default_display_colour_list();

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_LIST_HPP
