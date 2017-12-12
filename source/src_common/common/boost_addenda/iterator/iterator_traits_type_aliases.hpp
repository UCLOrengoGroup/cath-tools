/// \file
/// \brief The Boost Iterator traits type aliases header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_ITERATOR_ITERATOR_TRAITS_TYPE_ALIASES_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_ITERATOR_ITERATOR_TRAITS_TYPE_ALIASES_H

#include <boost/iterator/iterator_traits.hpp>

namespace cath {
	namespace common {

		template <typename T>
		using iterator_value_t      = typename boost::iterator_value     <T>::type;

		template <typename T>
		using iterator_reference_t  = typename boost::iterator_reference <T>::type;

		template <typename T>
		using iterator_pointer_t    = typename boost::iterator_pointer   <T>::type;

		template <typename T>
		using iterator_difference_t = typename boost::iterator_difference<T>::type;

		template <typename T>
		using iterator_category_t   = typename boost::iterator_category  <T>::type;

	} // namespace common
} // namespace cath

#endif
