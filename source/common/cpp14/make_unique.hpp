/// \file
/// \brief The make_unique header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_CPP14_MAKE_UNIQUE_H
#define _CATH_TOOLS_SOURCE_COMMON_CPP14_MAKE_UNIQUE_H

#include <memory>

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		template <typename T, typename... Args>
		std::unique_ptr<T> make_unique(Args&&... args
		                               ) {
			return std::unique_ptr<T>(
				new T( std::forward<Args>( args ) ... )
			);
		}
		
		/// \brief TODOCUMENT
		template <typename B, typename T, typename... Args>
		std::unique_ptr<B> make_base_unique(Args&&... args
		                                    ) {
			return std::unique_ptr<B>(
				new T( std::forward<Args>( args ) ... )
			);
		}
	} // namespace common
} // namespace cath

#endif
