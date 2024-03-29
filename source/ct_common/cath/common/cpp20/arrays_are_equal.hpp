/// \file
/// \brief The arrays_are_equal header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_ARRAYS_ARE_EQUAL_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_ARRAYS_ARE_EQUAL_HPP

#include <array>

namespace cath::common {

	template <typename T, size_t N>
	constexpr bool arrays_are_equal( const ::std::array<T, N> &prm_lhs, const ::std::array<T, N> &prm_rhs ) {
		for ( size_t ctr = 0; ctr < N; ++ctr ) {
			if ( prm_lhs[ ctr ] != prm_rhs[ ctr ] ) {
				return false;
			}
		}
		return true;
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_ARRAYS_ARE_EQUAL_HPP
