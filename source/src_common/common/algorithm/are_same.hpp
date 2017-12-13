/// \file
/// \brief The are_same() header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_ARE_SAME_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ALGORITHM_ARE_SAME_HPP

#include "common/cpp14/cbegin_cend.hpp"
#include "common/cpp17/invoke.hpp"
#include "common/function/ident.hpp"

namespace cath {
	namespace common {

		/// \brief Return whether all elements in the specified range are equal (evaluate true under the specified binary predicate)
		///        after being projected by the specified function, if any
		///
		///
		/// \todo Constrain this function so that it isn't ambiguous with the range version
		template <typename BeginItr,
		          typename EndItr,
		          typename Eq   = std::equal_to<>,
		          typename Proj = ident
		          >
		bool are_same_itr(BeginItr &&arg_begin_itr,     ///< The begin of the range to query
		                  EndItr   &&arg_end_itr,       ///< The end of the range to query
		                  Eq       &&arg_equal = Eq{},  ///< The binary predicate to use as the equality operator
		                  Proj     &&arg_proj  = Proj{} ///< The function to use to project the elements
		                  ) {
			if ( arg_begin_itr != arg_end_itr ) {
				const auto &first = invoke( arg_proj, *arg_begin_itr );
				for (auto start_itr = arg_begin_itr; start_itr != arg_end_itr; ++start_itr) {
					if ( ! invoke( arg_equal, first, invoke( arg_proj, *start_itr ) ) ) {
						return false;
					}
				}
			}
			return true;
		}

		/// \brief Return whether all elements in the specified range are equal (evaluate true under the specified binary predicate)
		///        after being projected by the specified function, if any
		template <typename Rng,
		          typename Eq   = std::equal_to<>,
		          typename Proj = ident
		          >
		bool are_same(Rng  &&arg_rng,           ///< The range to query
		              Eq   &&arg_equal = Eq{},  ///< The binary predicate to use as the equality operator
		              Proj &&arg_proj  = Proj{} ///< The function to use to project the elements
		              ) {
			return are_same_itr(
				common::cbegin      ( arg_rng   ),
				common::cend        ( arg_rng   ),
				std::forward< Eq   >( arg_equal ),
				std::forward< Proj >( arg_proj  )
			);
		}

	} // namespace common
} // namespace cath

#endif
