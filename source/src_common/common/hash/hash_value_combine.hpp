/// \file
/// \brief The hash_value_combine header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_HASH_HASH_VALUE_COMBINE_H
#define _CATH_TOOLS_SOURCE_COMMON_HASH_HASH_VALUE_COMBINE_H

#include <functional>

namespace cath {
	namespace common {

		/// \brief Add a new hash value into an existing hash seed
		inline void hash_value_combine(size_t       &arg_seed,      ///< The seed into which the hash value should be incorporated
		                               const size_t &arg_hash_value ///< The new hash value
		                               ) {
			arg_seed ^= arg_hash_value + 0x9e3779b9 + ( arg_seed << 6 ) + ( arg_seed >> 2 );
		}

	} // namespace common
} // namespace cath

#endif
