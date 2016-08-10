/// \file
/// \brief The constexpr_for_n class header

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

#ifndef CONSTEXPR_FOR_N_H_INCLUDED
#define CONSTEXPR_FOR_N_H_INCLUDED

#include "common/c++14/constexpr_min_max.h"

#include <utility>

namespace cath {
	namespace common {

		namespace detail {

			/// \brief Recursive implementation of constexpr_for_n to execute the function operator of F<I>()
			///        for all size_t values in half-open range [ 0, TO )
			template <template <size_t> class F, size_t FROM, size_t TO, typename ...Args>
			void constexpr_for_n_impl(Args && ...args) {
				// Execute the function operator of F<FROM>()
				F<FROM>()( std::forward<Args>( args )... );

				// If FROM + 1 is still less than TO, recursively call constexpr_for_n_impl() for FROM + 1
				if ( FROM + 1 < TO ) {
					constexpr_for_n_impl<F, cath::common::constexpr_min( FROM + 1, TO - 1 ), TO>( std::forward<Args>( args )... );
				}
			}

		}

		/// \brief Execute the function operator of F<I>() for all size_t values in half-open range [ 0, TO )
		///
		/// This is useful for executing some code for each entry of a std::array with the value as a constexpr
		///
		/// This is implemented in terms of the recursive constexpr_for_n_impl()
		template <template <size_t> class F, size_t TO,typename ...Args>
		void constexpr_for_n(Args && ...args) {
			detail::constexpr_for_n_impl<F, 0, TO>( std::forward<Args>( args )... );
		}

	}

}

#endif
