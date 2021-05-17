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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_STABLE_SORT_PROJ_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_STABLE_SORT_PROJ_HPP

#include <functional>

#include <boost/range/algorithm/stable_sort.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/function/ident.hpp"

namespace cath {
	namespace common {

		/// \brief Wrap boost::range::stable_sort() and add support for a projection function, a la Eric Niebler's range-v3 library
		///
		/// \todo Roll this out to simplify code
		template <typename Rng,
		          typename Pred = std::less<>,
		          typename Proj = ident>
		inline auto stable_sort_proj(Rng  &&prm_range,          ///< The range to stable_sort
		                             Pred &&prm_pred  = Pred{}, ///< The less-than predicate function
		                             Proj &&prm_proj  = Proj{}  ///< The projection function
		                             ) {
			return boost::range::stable_sort(
				prm_range,
				[&] (const auto &x, const auto &y) {
					return ::std::invoke(
						std::forward<Pred>( prm_pred ),
						::std::invoke(
							std::forward< Proj >( prm_proj ),
							x
						),
						::std::invoke(
							std::forward< Proj >( prm_proj ),
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
		inline auto stable_sort_proj_copy(Rng    prm_range,          ///< The range to stable_sort
		                                  Pred &&prm_pred  = Pred{}, ///< The less-than predicate function
		                                  Proj &&prm_proj  = Proj{}  ///< The projection function
		                                  ) {
			stable_sort_proj(
				prm_range,
				std::forward<Pred>( prm_pred ),
				std::forward<Proj>( prm_proj )
			);
			return prm_range;
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_STABLE_SORT_PROJ_HPP
