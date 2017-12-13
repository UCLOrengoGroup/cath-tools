/// \file
/// \brief The pair_hash class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_HASH_PAIR_HASH_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_HASH_PAIR_HASH_HPP

#include <functional>

namespace cath {
	namespace common {

		/// \brief Hashing for pairs of elements
		struct pair_hash {
		public:
			/// \brief Generic function operator for hashing a pair by
			///        combining the results of a standard hash of each value
			template <typename T, typename U>
			size_t operator()(const std::pair<T, U> &arg_pair ///< The
			                  ) const {
				return std::hash<T>()( arg_pair.first ) ^ std::hash<U>()( arg_pair.second );
			}
		};

	} // namespace common
} // namespace cath

#endif
