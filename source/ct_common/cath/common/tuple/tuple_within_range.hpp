/// \file
/// \brief The tuple_within_range header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_WITHIN_RANGE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_WITHIN_RANGE_HPP

#include <tuple>
#include <cstddef>

#include "cath/common/algorithm/variadic_and.hpp"
#include "cath/common/detail/make_static_const.hpp"
#include "cath/common/detail/tuple_index_sequence.hpp"

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Implementation function for tuple_within_range
			template <typename Tpl, size_t... Index>
			constexpr bool in_range_impl(const Tpl &prm_indexes,       ///< The index tuple whose values are to be checked
			                             const Tpl &prm_nums_cells,    ///< The tuple containing the values that the corresponding index values should be strictly less than
			                             std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                             ) {
				return variadic_and(
					std::get<Index>( prm_indexes ) >= 0 && std::get<Index>( prm_indexes ) < std::get<Index>( prm_nums_cells )...
				);
			}

			/// \brief Function object to check that each value in the specified index tuple is non-negative
			///        and is strictly less than the corresponding value in the specified nums tuple
			struct tuple_within_range_fn final {

				/// \brief Check that each value in the specified index tuple is non-negative
				///        and is strictly less than the corresponding value in the specified nums tuple
				template <typename Tpl>
				constexpr bool operator()(const Tpl &prm_indexes,   ///< The index tuple whose values are to be checked
				                          const Tpl &prm_nums_cells ///< The tuple containing the values that the corresponding index values should be strictly less than
				                          ) const {
					return in_range_impl(
						prm_indexes,
						prm_nums_cells,
						tuple_index_sequence<Tpl>{}
					);
				}

				tuple_within_range_fn()                              = delete;
				tuple_within_range_fn(const tuple_within_range_fn &) = delete;
				void operator=(const tuple_within_range_fn &)        = delete;
			};

		} // namespace detail

		MAKE_STATIC_CONST( detail::tuple_within_range_fn, tuple_within_range )

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_TUPLE_TUPLE_WITHIN_RANGE_HPP
