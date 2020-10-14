/// \file
/// \brief The layout_values_with_descs header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_LAYOUT_VALUES_WITH_DESCS_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_LAYOUT_VALUES_WITH_DESCS_HPP

#include <boost/optional.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/max_proj_element.hpp"
#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/type_aliases.hpp"

namespace cath {
	namespace common {

		/// \brief Return strings, one corresponding to each of the values in the specified range, in which the strings returned from
		///        the specified pair of functions for that value are laid out in columns with the specified separator
		template <typename Rng,
		          typename FnLhs,
		          typename FnRhs>
		str_vec layout_values_with_descs(const Rng         &prm_range,   ///< The range of values 
		                                 FnLhs              prm_fn_lhs,  ///< The first,  left-hand  function which must return a string when given a value of the range
		                                 FnRhs              prm_fn_rhs,  ///< The second, right-hand function which must return a string when given a value of the range
		                                 const std::string &prm_pair_sep ///< The separator with which to join the two functions' strings for each value
		                                 ) {
			const auto length_lhs_fn  = [&] (const auto &x) { return prm_fn_lhs( x ).length(); };
			const auto max_length_lhs = common::max_proj( prm_range, std::less<>{}, length_lhs_fn );

			return transform_build<str_vec>(
				prm_range,
				[&] (const auto &x) {
					return
						  prm_fn_lhs( x )
						+ std::string( max_length_lhs - length_lhs_fn( x ), ' ' ) + prm_pair_sep
						+ prm_fn_rhs( x );
				}
			);
		}

	} // namespace common
} // namespace cath

#endif
