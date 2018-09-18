/// \file
/// \brief The ident header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FUNCTION_IDENT_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FUNCTION_IDENT_HPP

#include <utility>

namespace cath {
	namespace common {

		/// \brief Function object to return the argument
		class ident final {
		public:
			/// \brief Return the argument using perfect-forwarding
			template <typename T>
			T && operator()(T &&prm_value
			                ) const noexcept {
				return std::forward<T>( prm_value );
			}

			/// \brief Indicate that this can be used in any smart, heterogeneous-lookup magic
			using is_transparent = void;
		};

	} // namespace common
} // namespace cath

#endif
