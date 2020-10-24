/// \file
/// \brief The display_colour_gradient class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_GRADIENT_HPP
#define _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_GRADIENT_HPP

#include "cath/display_colour/display_colour.hpp"
#include "cath/display_colour/display_colour_type_aliases.hpp"

#include <cstddef>
#include <vector>

namespace cath {

	/// \brief TODOCUMENT
	class display_colour_gradient final {
	private:
		/// \brief TODOCUMENT
		display_colour_vec colour_points;

		/// \brief TODOCUMENT
		size_t             steps_in_between_points;

		void check_values() const;

	public:
		display_colour_gradient(display_colour_vec,
		                        const size_t &);

		size_t get_steps_in_between_points() const;
		size_t get_num_colour_points() const;
		const display_colour & get_colour_point_of_index(const size_t &) const;
	};

	size_t get_num_colours(const display_colour_gradient &);
	display_colour get_colour_of_fraction(const display_colour_gradient &,
	                                      const double &);
	display_colour get_colour_of_index(const display_colour_gradient &,
	                                   const size_t &);
	display_colour_gradient make_default_light_colour_gradient();
	display_colour_gradient make_default_colour_gradient();
	display_colour_gradient make_default_dark_colour_gradient();
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_GRADIENT_HPP
