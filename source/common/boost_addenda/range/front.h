/// \file
/// \brief The front header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_FRONT_H
#define _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_FRONT_H

#include <boost/range.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.h"
#include "common/cpp14/cbegin_cend.h"

namespace cath {
	namespace common {

		/// \brief Return a non-const reference to the first element of a range
		///
		/// This is a non-member function for ranges that don't provide their own front method.
		/// This is useful for non-member ranges
		template <typename T>
		inline range_reference_t<T> front(T &arg_range ///< The range to query
		                                  ) {
			return * std::begin( arg_range );
		}

		/// \brief Return a const reference to the first element of a range
		///
		/// \copydetails front()
		template <typename T>
		inline range_const_reference_t<T> front(const T &arg_range ///< The range to query
		                                        ) {
			return * common::cbegin( arg_range );
		}

	} // namespace common
} // namespace cath

#endif
