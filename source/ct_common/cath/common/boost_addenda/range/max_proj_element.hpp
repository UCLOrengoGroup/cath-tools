/// \file
/// \brief The max_proj_element header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_MAX_PROJ_ELEMENT_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_MAX_PROJ_ELEMENT_HPP

#include <boost/range/algorithm/max_element.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/function/ident.hpp"

namespace cath {
	namespace common {

		/// \brief Wrap boost::range::max_element() and add support for a projection function, a la Eric Niebler's range-v3 library
		///
		/// \todo Roll this out to simplify code
		template <typename Rng,
		          typename Pred = std::less<>,
		          typename Proj = ident>
		inline auto max_proj_element(Rng  &&prm_range,          ///< The range to query
		                             Pred &&prm_pred  = Pred{}, ///< The less-than predicate function
		                             Proj &&prm_proj  = Proj{}  ///< The projection function
		                             ) {
			return boost::range::max_element(
				prm_range,
				[&] (const auto &x, const auto &y) {
					return std::forward<Pred>( prm_pred )( prm_proj( x ), prm_proj( y ) );
				}
			);
		}

		/// \brief Return the projection of the max_proj_element() of the specified range
		///
		/// \todo Roll this out to simplify code
		template <typename Rng,
		          typename Pred = std::less<>,
		          typename Proj = ident>
		inline auto max_proj(Rng  &&prm_range,          ///< The range to query
		                     Pred &&prm_pred  = Pred{}, ///< The less-than predicate function
		                     Proj &&prm_proj  = Proj{}  ///< The projection function
		                     ) {
			return prm_proj( *max_proj_element(
				std::forward< Rng  >( prm_range ),
				std::forward< Pred >( prm_pred  ),
				std::forward< Proj >( prm_proj  )
			) );
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_MAX_PROJ_ELEMENT_HPP
