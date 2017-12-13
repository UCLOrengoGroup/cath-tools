/// \file
/// \brief The debug_numeric_cast header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_DEBUG_NUMERIC_CAST_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_DEBUG_NUMERIC_CAST_HPP

#include <boost/numeric/conversion/cast.hpp>

namespace cath {

	/// \brief TODOCUMENT
	template <typename Target, typename Source>
	inline Target debug_unwarned_numeric_cast(const Source &arg_value ///< TODOCUMENT
	                                          ) {
#ifndef NDEBUG
		return boost::numeric_cast<Target>( arg_value );
#else
		return static_cast<Target>( arg_value );
#endif
	}

	/// \brief TODOCUMENT
	template <typename Target, typename Source>
	inline Target debug_numeric_cast(const Source &arg_value ///< TODOCUMENT
	                                 ) {
		using stripped_target  = typename std::remove_reference<typename std::remove_const<Target>::type>::type;
		using stripped_source  = typename std::remove_reference<typename std::remove_const<Source>::type>::type;
		static_assert(
			! std::is_same<stripped_target, stripped_source>::value,
			"debug_numeric_cast is being used to cast a type to itself"
		);
		return debug_unwarned_numeric_cast<Target, Source>( arg_value );
	}
} // namespace cath

#endif
