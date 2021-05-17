/// \matrix
/// \brief The matrix_index header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_MATRIX_MATRIX_INDEX_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_MATRIX_MATRIX_INDEX_HPP

#include <type_traits>

namespace cath {
	namespace common {

		/// \brief Get the (0-offset) index of the specified element in the strict upper half of the square matrix
		///        of specified size (where the counting is along the 0-th row, then the 1-st row etc)
		///
		/// For instance, the following will print out the numbers 0 through 9 :
		///
		/// \code{.cpp}
		///	constexpr int N = 5;
		/// for (const int &i : indices( N ) ) {
		///   for (const int &i : irange( i + 1, N ) ) {
		///     std::cerr << get_zero_index_of_strict_upper_half_matrix( i, j, N ) << "\n";
		///   }
		/// }
		/// \endcode
		///
		/// ...which corresponds to this:
		///
		/// \verbatim
		/// .  0  1  2  3
		/// .  .  4  5  6
		/// .  .  .  7  8
		/// .  .  .  .  9
		/// .  .  .  .  .
		/// \endverbatim
		///
		/// The formula being used is: \f$ j + \frac{ i . \left( 2n -3 -i \right) }{2} \f$
		template <typename Int>
		constexpr Int get_zero_index_of_strict_upper_half_matrix(const Int &prm_row, ///< The row index, in \f$ [ 0,   	    arg\_n -1 ) \f$
		                                                         const Int &prm_col, ///< The col_index, in \f$ ( arg\_row, arg\_n    ) \f$
		                                                         const Int &prm_n    ///< The width/height of the square matrix
		                                                         ) {
			static_assert( std::is_integral_v<Int>, "The type used for matrix calculations must be an integral type" );
			return (
				prm_col
				+
				(
					prm_row
					*
					(
						( static_cast<Int>( 2 ) * prm_n )
						-
						static_cast<Int>( 3 )
						-
						prm_row
					)
					/ static_cast<Int>( 2 )
				)
				-
				1
			);
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_MATRIX_MATRIX_INDEX_HPP
