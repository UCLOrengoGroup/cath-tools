/// \file
/// \brief The minmax_element header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_MINMAX_ELEMENT_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_MINMAX_ELEMENT_HPP

#include <boost/algorithm/minmax_element.hpp>

#include "cath/common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		template <typename RNG, typename FN>
		inline auto minmax_element(RNG &prm_range ///< TODOCUMENT
		                           ) {
			return boost::minmax_element(
				common::cbegin( prm_range ),
				common::cend  ( prm_range )
			);
		}

		/// \brief TODOCUMENT
		template <typename RNG, typename FN>
		inline auto minmax_element(RNG &prm_range,   ///< TODOCUMENT
		                           FN   prm_function ///< TODOCUMENT
		                           ) {
			return boost::minmax_element(
				common::cbegin( prm_range ),
				common::cend  ( prm_range ),
				prm_function
			);
		}

	} // namespace common
} // namespace cath

#endif
