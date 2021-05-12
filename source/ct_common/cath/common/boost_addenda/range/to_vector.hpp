/// \file
/// \brief The to_vector header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_TO_VECTOR_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_TO_VECTOR_HPP

#include <iterator>
#include <string>
#include <vector>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/detail/make_static_const.hpp"

namespace cath::common {

	namespace detail {

		template <template <typename...> typename Cont>
		struct to_container_of_fn {
			template <typename Rng>
			Cont<range_value_t<Rng>> operator()( Rng &&prm_rng ) const {
				return Cont<range_value_t<Rng>>( ::std::cbegin( prm_rng ), ::std::cend( prm_rng ) );
			}
		};

		template <template <typename...> typename Cont, typename Rng>
		auto operator|( Rng &&prm_rng, const to_container_of_fn<Cont> &prm_to_container_of_fn ) {
			return prm_to_container_of_fn( ::std::forward<Rng>( prm_rng ) );
		}

	} // namespace detail

	/// Make a string from a range
	///
	/// Use by piping the range through to_string
	constexpr detail::to_container_of_fn<::std::basic_string> to_string;

	/// Make a vector from a range
	///
	/// Use by piping the range through to_string
	constexpr detail::to_container_of_fn<::std::vector> to_vector;

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_TO_VECTOR_HPP
