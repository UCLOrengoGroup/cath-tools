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

#ifndef _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_HPP
#define _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_HPP

#include <cmath>
#include <string>
#include <vector>

#include <boost/operators.hpp>
#include <boost/throw_exception.hpp>

#include <fmt/core.h>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/make_type_of_first_n.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {

	/// \brief TODOCUMENT
	///
	/// \todo Consider making this constexpr
	class display_colour final : private boost::equivalent<display_colour, boost::totally_ordered<display_colour>> {
		/// \brief TODOCUMENT
		double r = 0.0;

		/// \brief TODOCUMENT
		double g = 0.0;

		/// \brief TODOCUMENT
		double b = 0.0;

		static constexpr void check_component_value( const double & );

	  public:
		/// \brief Default ctor for display_colour (makes BLACK)
		constexpr display_colour() = default;
		constexpr display_colour( const double &, const double &, const double & );

		[[nodiscard]] constexpr const double &get_r() const;
		[[nodiscard]] constexpr const double &get_g() const;
		[[nodiscard]] constexpr const double &get_b() const;

		/// \brief TODOCUMENT
		static constexpr ::std::string_view COMPONENT_SEPARATOR{ "," };

		/// \brief Return whether the first colour is less than the second
		///
		/// This is a simple ordering that compares on the red, green and then blue components
		///
		/// \todo Is this actually useful? If so, document the context; if not, consider removing.
		///
		/// \relates display_colour
		friend constexpr bool operator<( const display_colour &prm_lhs, ///< The first  display_colour to compare
		                                 const display_colour &prm_rhs  ///< The second display_colour to compare
		                                 ) {
			return ( ::std::tie( prm_lhs.get_r(), prm_lhs.get_g(), prm_lhs.get_b() )
			         < ::std::tie( prm_rhs.get_r(), prm_rhs.get_g(), prm_rhs.get_b() ) );
		}
	};

	/// \brief TODOCUMENT
	constexpr void display_colour::check_component_value( const double &prm_component_value ///< TODOCUMENT
	                                                      ) {
		if ( !( prm_component_value >= 0.0 && prm_component_value <= 1.0 ) ) {
			BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
			  ::fmt::format( "Viewer colour component value {} is not between 0 and 1", prm_component_value ) ) );
		}
	}

	/// \brief Ctor for display_colour
	constexpr display_colour::display_colour( const double &prm_r, ///< TODOCUMENT
	                                          const double &prm_g, ///< TODOCUMENT
	                                          const double &prm_b  ///< TODOCUMENT
	                                          ) :
	        r( prm_r ), g( prm_g ), b( prm_b ) {
		check_component_value( r );
		check_component_value( g );
		check_component_value( b );
	}

	/// \brief Getter for the red component (between 0 and 1)
	constexpr const double &display_colour::get_r() const {
		return r;
	}

	/// \brief Getter for the green component (between 0 and 1)
	constexpr const double &display_colour::get_g() const {
		return g;
	}

	/// \brief Getter for the blue component (between 0 and 1)
	constexpr const double &display_colour::get_b() const {
		return b;
	}

	/// \brief TODOCUMENT
	///
	/// \relates display_colour
	///
	/// \param prm_components TODOCUMENT
	constexpr display_colour display_colour_from_components( const ::std::array<double, 3> &prm_components ) {
		return common::make_type_of_first_n<display_colour, 3>( prm_components );
	}

	display_colour display_colour_from_string( const std::string & );
	std::string    hex_string_of_colour( const display_colour & );
	std::string    comma_separated_string_of_display_colour( const display_colour & );
	std::ostream & operator<<( std::ostream &, const display_colour & );

	/// \brief TODOCUMENT
	inline constexpr display_colour BLACK{ 0.00, 0.00, 0.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour WHITE        { 1.00, 1.00, 1.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour RED          { 1.00, 0.00, 0.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour LIGHT_RED    { 1.00, 0.30, 0.30 };

	/// \brief TODOCUMENT
	inline constexpr display_colour YELLOW       { 1.00, 1.00, 0.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour DARK_YELLOW  { 0.70, 0.70, 0.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour GREEN        { 0.00, 1.00, 0.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour LIGHT_GREEN  { 0.30, 1.00, 0.30 };

	/// \brief TODOCUMENT
	inline constexpr display_colour CYAN         { 0.00, 1.00, 1.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour DARK_CYAN    { 0.00, 0.70, 0.70 };

	/// \brief TODOCUMENT
	inline constexpr display_colour BLUE         { 0.00, 0.00, 1.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour LIGHT_BLUE   { 0.30, 0.30, 1.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour MAGENTA      { 1.00, 0.00, 1.00 };

	/// \brief TODOCUMENT
	inline constexpr display_colour DARK_MAGENTA { 0.70, 0.00, 0.70 };

	/// \brief TODOCUMENT
	///
	/// \TODO: Make constexpr
	///
	/// \relates display_colour
	constexpr display_colour rgb_mid_point( const display_colour &prm_colour1,      ///< TODOCUMENT
	                                        const display_colour &prm_colour2,      ///< TODOCUMENT
	                                        const double &        prm_fraction_thru ///< TODOCUMENT
	                                        ) {
		if ( !( prm_fraction_thru >= 0.0 && prm_fraction_thru <= 1.0 ) ) {
			BOOST_THROW_EXCEPTION(
			  common::invalid_argument_exception( "Fraction through must be a finite number between 0 and 1" ) );
		}

		const double fraction_from = 1.0 - prm_fraction_thru;
		return display_colour( fraction_from * prm_colour1.get_r() + prm_fraction_thru * prm_colour2.get_r(),
		                       fraction_from * prm_colour1.get_g() + prm_fraction_thru * prm_colour2.get_g(),
		                       fraction_from * prm_colour1.get_b() + prm_fraction_thru * prm_colour2.get_b() );
	}

	/// \brief Return a copy of the specified colour that has been lightened by the specified fraction
	///
	/// \relates display_colour
	constexpr display_colour lighten_by_fraction( const display_colour &prm_colour, ///< The colour from which to start
	                                              const double &prm_fraction ///< The fraction change (0.0 gives the original colour; 1.0 gives white)
	                                              ) {
		return rgb_mid_point( prm_colour, WHITE, prm_fraction );
	}

	/// \brief Return a copy of the specified colour that has been darkened by the specified fraction
	///
	/// \relates display_colour
	constexpr display_colour darken_by_fraction( const display_colour &prm_colour, ///< The colour from which to start
	                                             const double &prm_fraction ///< The fraction change (0.0 gives the original colour; 1.0 gives black)
	                                             ) {
		return rgb_mid_point( prm_colour, BLACK, prm_fraction );
	}

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_DISPLAY_COLOUR_CATH_DISPLAY_COLOUR_DISPLAY_COLOUR_HPP
