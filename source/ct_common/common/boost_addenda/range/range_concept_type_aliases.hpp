/// \file
/// \brief The Boost Range concept type aliases header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_RANGE_CONCEPT_TYPE_ALIASES_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_RANGE_RANGE_CONCEPT_TYPE_ALIASES_HPP

#include <boost/range/metafunctions.hpp>

#include "common/cpp17/void_t.hpp"

#include <vector>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Primary template for implementation of range_const_reference
			///
			/// By default, use range_reference<const T>
			template <typename T, typename = detail::void_t<> >
			struct range_const_reference_impl {
				using type = typename boost::range_reference<const T>::type;
			};

			/// \brief Specialization of template for implementation of range_const_reference
			///
			/// If T has a const_reference member type, then prefer to use that
			///
			/// This uses the void_t trick to remove this specialisation if T doesn't contain a const_reference member type.
			template <typename T>
			struct range_const_reference_impl< T, detail::void_t< typename T::const_reference > > {
				using type = typename T::const_reference;
			};

		} // namespace detail

		template <typename T>
		using range_iterator_t               = typename boost::range_iterator<T>::type;

		template <typename T>
		using range_const_iterator_t         = typename boost::range_iterator<const T>::type;

		template <typename T>
		using range_value_t                  = typename boost::range_value<T>::type;

		template <typename T>
		using range_reference_t              = typename boost::range_reference<T>::type;

		template <typename T>
		using range_const_reference_t        = typename detail::range_const_reference_impl<T>::type;

		template <typename T>
		using range_pointer_t                = typename boost::range_pointer<T>::type;

		template <typename T>
		using range_category_t               = typename boost::range_category<T>::type;

		template <typename T>
		using range_difference_t             = typename boost::range_difference<T>::type;

		template <typename T>
		using range_reverse_iterator_t       = typename boost::range_reverse_iterator<T>::type;

		template <typename T>
		using range_reverse_const_iterator_t = typename boost::range_reverse_iterator<const T>::type;

		template <typename T>
		using has_range_iterator_t           = typename boost::has_range_iterator<T>::type;

		template <typename T>
		using has_range_const_iterator_t     = typename boost::has_range_const_iterator<T>::type;

	} // namespace common
} // namespace cath

#endif
