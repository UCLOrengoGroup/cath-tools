/// \file
/// \brief The cbegin() / crbegin() / cend() / crend() header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CPP14_CBEGIN_CEND_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CPP14_CBEGIN_CEND_HPP

#include <iterator>

namespace cath {
	namespace common {

		/// \brief Temporary replacement for C++14's std::cbegin()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto cbegin(const T &arg_range
		                   )->decltype( std::begin( arg_range ) ) {
			return std::begin( arg_range );
		}

		/// \brief Temporary replacement for C++14's std::cend()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto cend(const T &arg_range
		                 )->decltype( std::end( arg_range ) ) {
			return std::end( arg_range );
		}

		/// \brief Temporary replacement for C++14's non-const overload of std::rbegin()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto rbegin(T& arg_range
		                   ) -> decltype( arg_range.rbegin() ) {
			return arg_range.rbegin();
		}

		/// \brief Temporary replacement for C++14's const overload of std::rbegin()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto rbegin(const T& arg_range
		                   ) -> decltype( arg_range.rbegin() ) {
			return arg_range.rbegin();
		}

		/// \brief Temporary replacement for C++14's non-const overload of std::rend()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto rend(T& arg_range
		                 ) -> decltype( arg_range.rend() ) {
			return arg_range.rend();
		}

		/// \brief Temporary replacement for C++14's const overload of std::rend()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto rend(const T &arg_range
		                 ) -> decltype( arg_range.rend() ) {
			return arg_range.rend();
		}


		/// \brief Temporary replacement for C++14's const overload of std::crbegin()
		///
		/// \todo Come C++14, remove this and use the proper version instead
		template <typename T>
		inline auto crbegin(const T &arg_range
		                    ) -> decltype( ::cath::common::rbegin( arg_range ) ) {
			return ::cath::common::rbegin( arg_range );
		}

		template <typename T>
		inline auto crend(const T& arg_range
		                  ) -> decltype( ::cath::common::rend( arg_range ) ) {
			return ::cath::common::rend( arg_range );
		}

	} // namespace common
} // namespace cath

#endif
