/// \file
/// \brief The sort_proj header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_SORT_PROJ_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_SORT_PROJ_HPP

#include <boost/range/algorithm/sort.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/cpp17/invoke.hpp"
#include "common/function/ident.hpp"

namespace cath {
	namespace common {

		/// \brief Wrap boost::range::sort() and add support for a projection function, a la Eric Niebler's range-v3 library
		///
		/// \todo Roll this out to simplify code
		template <typename Rng,
		          typename Pred = std::less<>,
		          typename Proj = ident>
		inline auto sort_proj(Rng  &&arg_range,          ///< The range to sort
		                      Pred &&arg_pred  = Pred{}, ///< The less-than predicate function
		                      Proj &&arg_proj  = Proj{}  ///< The projection function
		                      ) {
			return boost::range::sort(
				arg_range,
				[&] (const auto & x, const auto & y) {
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

	} // namespace common
} // namespace cath

#endif
