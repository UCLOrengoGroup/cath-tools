/// \file
/// \brief The sort_copy header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#ifndef SORT_COPY_H_INCLUDED
#define SORT_COPY_H_INCLUDED

#include <boost/range/algorithm/sort.hpp>

namespace cath {
	namespace common {

		/// \brief Convenience function for making a sorted copy of a range
		template <typename R>
		R sort_copy(R arg_range ///< The range on which the sorted copy should be based
		            ) {
			boost::range::sort( arg_range );
			return arg_range;
		}

		/// \overload
		template <typename R, typename P>
		R sort_copy(R arg_range,   ///< The range on which the sorted copy should be based
		            P arg_bin_pred ///< The binary predicate to use as a less-than operator for sorting
		            ) {
			boost::range::sort( arg_range, arg_bin_pred );
			return arg_range;
		}

	}
}

#endif
