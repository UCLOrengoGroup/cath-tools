/// \file
/// \brief The cath_to_string header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_CATH_TO_STRING_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_CATH_TO_STRING_HPP

#include <array>
#include <optional>
#include <string>
#include <tuple>

#include <boost/algorithm/string/join.hpp>
#include <boost/core/demangle.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <fmt/core.h>

#include "cath/common/type_traits.hpp"
#include "cath/common/type_traits/is_template_of_type.hpp"
#include "cath/common/type_traits/is_tuple.hpp"

namespace cath::common {

	template <typename T>
	::std::string cath_to_string( T &&prm_value ) {
		if constexpr ( is_same_modulo_cvref_v<T, ::std::string> ) {
			return ::std::forward<T>( prm_value );
		} else if constexpr ( is_same_modulo_cvref_v<T, bool> ) {
			return prm_value ? "true" : "false";
		} else if constexpr ( is_same_modulo_cvref_v<T, ::std::nullopt_t> ) {
			return "nullopt";
		} else if constexpr ( is_template_of_type_v<common::remove_cvref_t<T>, ::std::optional> ) {
			return ::fmt::format( "std::optional<{}>({})",
			                      ::boost::core::demangle( typeid( typename common::remove_cvref_t<T>::value_type ).name() ),
			                      prm_value.has_value() ? cath_to_string( *prm_value ) : "nullopt" );
		} else if constexpr ( is_tuple_mod_cvref_v<T> || is_template_of_type_v<common::remove_cvref_t<T>, ::std::pair> ) {
			// clang-format off
			return ::std::apply(
				[]( const auto &...vals ) {
					return ::fmt::format(
						"std::{}<{}>({})",
						is_tuple_mod_cvref_v<T> ? "tuple" : "pair",
						::boost::algorithm::join(
							::std::array{ ::boost::core::demangle( typeid( decltype( vals ) ).name() )... },
							", "
						),
						::boost::algorithm::join(
							::std::array{ cath_to_string( vals )...  },
							", "
						)
					);
				},
				prm_value
			);
			// clang-format on
		} else {
			using ::std::to_string;
			return to_string( ::std::forward<T>( prm_value ) );
		}
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_STRING_CATH_TO_STRING_HPP
