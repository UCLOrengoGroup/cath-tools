/// \file
/// \brief The temporary check_offset_1 header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TEMP_CHECK_OFFSET_1_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TEMP_CHECK_OFFSET_1_HPP

#include "cath/common/config.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

#include <cstddef>

namespace cath {
	/// \brief A helper function to check offset_1 values are not 0
	///
	/// This should be temporary whilst transitioning from some old offset_1 code
	inline void check_offset_1(const size_t &prm_index_offset_1 ///< The index (offset 1)
	                           ) {
		if constexpr ( common::IS_IN_DEBUG_MODE ) {
			if ( prm_index_offset_1 == 0 ) {
				BOOST_THROW_EXCEPTION(
				  cath::common::invalid_argument_exception( "Index specified with offset of 1 cannot be 0" ) );
			}
		}
	}
} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TEMP_CHECK_OFFSET_1_HPP
