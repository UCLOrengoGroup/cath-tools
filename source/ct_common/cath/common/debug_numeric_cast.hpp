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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_DEBUG_NUMERIC_CAST_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_DEBUG_NUMERIC_CAST_HPP

#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/config.hpp"
#include "cath/common/type_traits.hpp"

namespace cath {

	/// \brief TODOCUMENT
	template <typename Target, typename Source>
	constexpr Target debug_unwarned_numeric_cast(const Source &prm_value ///< TODOCUMENT
	                                             ) {
		if constexpr ( common::IS_IN_DEBUG_MODE ) {
			return ::boost::numeric_cast<Target>( prm_value );
		} else {
			return static_cast<Target>( prm_value );
		}
	}

	/// \brief TODOCUMENT
	template <typename Target, typename Source>
	constexpr Target debug_numeric_cast(const Source &prm_value ///< TODOCUMENT
	                                    ) {
		static_assert(
			! common::is_same_modulo_cvref_v<Target, Source>,
			"debug_numeric_cast is being used to cast a type to itself"
		);
		return debug_unwarned_numeric_cast<Target, Source>( prm_value );
	}

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_DEBUG_NUMERIC_CAST_HPP
