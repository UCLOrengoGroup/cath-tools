/// \file
/// \brief The more header

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

#ifndef MORE_H_INCLUDED
#define MORE_H_INCLUDED

#include <utility>

namespace cath {
	namespace common {

		/// \brief Function object to return the result of the less-than operation in the opposite direction
		class more final {
		public:
			/// \brief Return the argument using perfect-forwarding
			template <typename T, typename U>
			constexpr bool operator()(T &&arg_lhs, ///< TODOCUMENT
			                          U &&arg_rhs  ///< TODOCUMENT
			                          ) const {
				return std::forward<U>( arg_rhs ) < std::forward<T>( arg_lhs );
			}

			/// \brief Indicate that this can be used in any smart, heterogeneous-lookup magic
			using is_transparent = void;
		};

	}
}

#endif
