/// \file
/// \brief The lexical_cast_itr class header

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

#ifndef LEXICAL_CAST_ITR_H_INCLUDED
#define LEXICAL_CAST_ITR_H_INCLUDED

#include <boost/iterator/transform_iterator.hpp>
#include <boost/lexical_cast.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.h"

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Simple functor to lexical_cast from type F to type T
			template <typename F, typename T>
			class lexical_cast_value final {
			public:
				T operator()(const F &) const;
			};

			/// \brief Function operator to lexical_cast from type F to type F
			template <typename F, typename T>
			T lexical_cast_value<F, T>::operator()(const F &arg_value ///< The value to lexical_cast
			                                       ) const {
				return boost::lexical_cast<T>( arg_value );
			}
		}

		/// \brief Type alias for a lexical_cast_itr for a given destination type, T, and range, RNG.
		///
		/// For a const_iterator, use a const RNG type.
		///
		/// This  is implemented as a boost::transform_iterator<> using lexical_cast_value
		template <typename T, typename RNG>
		using lexical_cast_itr = boost::transform_iterator<
			detail::lexical_cast_value<range_value_t<RNG>, T>,
			range_iterator_t<RNG>
		>;

	}
}

#endif