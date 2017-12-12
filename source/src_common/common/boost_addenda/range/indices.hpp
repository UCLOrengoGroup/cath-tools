/// \file
/// \brief The indices header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_INDICES_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_INDICES_H

#include <boost/range/irange.hpp>

#include "common/boost_addenda/range/indices.hpp"

namespace cath {
	namespace common{

		/// \brief Return an integer_range between zero and the specified value
		template <typename T>
		::boost::integer_range<T> indices(const T &arg_n ///< The (one-past-the) end value
		                                  ) {
			return ::boost::irange( static_cast<T>( 0 ), arg_n );
		}

	}
}

#endif
