/// \file
/// \brief The constexpr_ignore_unused header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CONSTEXPR_IGNORE_UNUSED_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CONSTEXPR_IGNORE_UNUSED_HPP

namespace cath {
	namespace common {

		/// \brief Like Boost's ignore unused but constexpr
		///
		/// See http://www.boost.org/libs/core/doc/html/core/ignore_unused.html
		template <typename... Ts>
		constexpr bool constexpr_ignore_unused(const Ts &...) {
			return true;
		}
		
	} // namespace common
} // namespace cath

#endif
