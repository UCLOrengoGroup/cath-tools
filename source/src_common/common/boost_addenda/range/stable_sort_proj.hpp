/// \file
/// \brief The stable_sort_proj header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_STABLE_SORT_PROJ_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_STABLE_SORT_PROJ_H

#include <boost/range/algorithm/stable_sort.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/cpp17/invoke.hpp"
#include "common/function/ident.hpp"

namespace cath {
	namespace common {

		/// \brief Wrap boost::range::stable_sort() and add support for a projection function, a la Eric Niebler's range-v3 library
		///
		/// \todo Roll this out to simplify code
		template <typename Rng,
		          typename Pred = std::less<>,
		          typename Proj = ident>
		inline auto stable_sort_proj(Rng  &&arg_range,          ///< The range to stable_sort
		                             Pred &&arg_pred  = Pred{}, ///< The less-than predicate function
		                             Proj &&arg_proj  = Proj{}  ///< The projection function
		                             ) {
			/// \todo Come Clang (>= 3.7?) with fix, drop this type alias and use generic lambdas
			///
			/// This is a bit of a mess. Unlike sort(), stable_sort() appears to be unhappy if the
			/// predicate arguments are non-const lvalue references because it tries to pass const values.
			///
			/// (This seems to contradict the specification on cppreference - I hope to report that conflict.)
			///
			/// Unfortunately, range_const_reference_t<> doesn't return the correct type from Boost Range
			/// adaptors (even after I improved it to specifically detect and use const_reference member types).
			///
			/// So here, the following adds const to lvalue references (which requires stripping the lvalue reference
			/// off, adding the const and replacing the lvalue reference because you can't add const directly to the lvalue
			/// reference)
			using const_reference_type_naive = common::range_const_reference_t<Rng>;
			using const_reference_type       = std::conditional_t<
				std::is_lvalue_reference<const_reference_type_naive>::value,
				std::add_lvalue_reference_t<std::add_const_t<std::decay_t<const_reference_type_naive>>>,
				                            std::add_const_t<std::decay_t<const_reference_type_naive>>
			>;

			return boost::range::stable_sort(
				arg_range,
				[&] (const_reference_type x, const_reference_type y) {
					return common::invoke(
						std::forward<Pred>( arg_pred ),
						common::invoke(
							std::forward< Proj >( arg_proj ),
							x
						),
						common::invoke(
							std::forward< Proj >( arg_proj ),
							y
						)
					);
				}
			);
		}

		/// \brief Return a copy of the specified range after stably-sorting using the specified less-than predicate and projection function
		///
		/// \todo Roll this out to simplify code
		template <typename Rng,
		          typename Pred = std::less<>,
		          typename Proj = ident>
		inline auto stable_sort_proj_copy(Rng    arg_range,          ///< The range to stable_sort
		                                  Pred &&arg_pred  = Pred{}, ///< The less-than predicate function
		                                  Proj &&arg_proj  = Proj{}  ///< The projection function
		                                  ) {
			stable_sort_proj(
				arg_range,
				std::forward<Pred>( arg_pred ),
				std::forward<Proj>( arg_proj )
			);
			return arg_range;
		}

	} // namespace common
} // namespace cath

#endif
