/// \file
/// \brief The sort_copy header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_COPY_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_COPY_HPP

#include <boost/range/algorithm/sort.hpp>

namespace cath {
	namespace common {

		/// \brief Convenience function for making a sorted copy of a range
		template <typename R>
		R sort_copy(R prm_range ///< The range on which the sorted copy should be based
		            ) {
			boost::range::sort( prm_range );
			return prm_range;
		}

		/// \overload
		template <typename R, typename P>
		R sort_copy(R prm_range,   ///< The range on which the sorted copy should be based
		            P prm_bin_pred ///< The binary predicate to use as a less-than operator for sorting
		            ) {
			boost::range::sort( prm_range, prm_bin_pred );
			return prm_range;
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_ALGORITHM_SORT_COPY_HPP
