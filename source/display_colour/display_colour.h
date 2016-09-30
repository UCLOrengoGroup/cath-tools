/// \file
/// \brief The display_colour class header

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

#ifndef DISPLAY_COLOUR_H_INCLUDED
#define DISPLAY_COLOUR_H_INCLUDED

#include <boost/operators.hpp>

#include "common/type_aliases.h"

#include <string>
#include <vector>

namespace cath {

	/// \brief TODOCUMENT
	///
	/// \todo Consider making this constexpr
	class display_colour final : private boost::equivalent      < display_colour,
	                                     boost::totally_ordered < display_colour > > {
		/// \brief TODOCUMENT
		double r = 0.0;

		/// \brief TODOCUMENT
		double g = 0.0;

		/// \brief TODOCUMENT
		double b = 0.0;

		static void check_component_value(const double &);

	public:
		/// \brief Default ctor for display_colour (makes BLACK)
		display_colour() = default;
		display_colour(const double &,
		               const double &,
		               const double &);

		const double & get_r() const;
		const double & get_g() const;
		const double & get_b() const;

		static const std::string    COMPONENT_SEPARATOR;
		static const display_colour BLACK;
		static const display_colour WHITE;

		static const display_colour RED;
		static const display_colour LIGHT_RED;
		static const display_colour YELLOW;
		static const display_colour DARK_YELLOW;
		static const display_colour GREEN;
		static const display_colour LIGHT_GREEN;
		static const display_colour CYAN;
		static const display_colour DARK_CYAN;
		static const display_colour BLUE;
		static const display_colour LIGHT_BLUE;
		static const display_colour MAGENTA;
		static const display_colour DARK_MAGENTA;
	};

	display_colour lighten_by_fraction(const display_colour &,
	                                   const double &);
	display_colour darken_by_fraction(const display_colour &,
	                                  const double &);

	bool operator<(const display_colour &,
	               const display_colour &);

	display_colour display_colour_from_string(const std::string &);
	display_colour display_colour_from_components(const doub_vec &);
	std::string    hex_string_of_colour(const display_colour &);
	std::string    comma_separated_string_of_display_colour(const display_colour &);
	display_colour rgb_mid_point(const display_colour &,
	                             const display_colour &,
	                             const double &);
	std::ostream & operator<<(std::ostream &,
	                          const display_colour &);
}

#endif
