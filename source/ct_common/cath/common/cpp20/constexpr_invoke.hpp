/// \file
/// \brief The constexpr_invoke header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_CONSTEXPR_INVOKE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_CONSTEXPR_INVOKE_HPP

#include <type_traits>
#include <utility>

// Implementation copied/pasted from https://en.cppreference.com/w/cpp/utility/functional/invoke
//                                                                                (17th May 2021)

namespace cath::common {
	namespace detail {

		template <class>
		constexpr bool is_reference_wrapper_v = false;
		template <class U>
		constexpr bool is_reference_wrapper_v<::std::reference_wrapper<U>> = true;

		template <class T, class Type, class T1, class... Args>
		constexpr decltype( auto ) constexpr_invoke_impl( Type T::*f, T1 &&t1, Args &&...args ) {
			if constexpr ( ::std::is_member_function_pointer_v<decltype( f )> ) {
				if constexpr ( ::std::is_base_of_v<T, ::std::decay_t<T1>> )
					return ( ::std::forward<T1>( t1 ).*f )( ::std::forward<Args>( args )... );
				else if constexpr ( is_reference_wrapper_v<::std::decay_t<T1>> )
					return ( t1.get().*f )( ::std::forward<Args>( args )... );
				else
					return ( ( *::std::forward<T1>( t1 ) ).*f )( ::std::forward<Args>( args )... );
			} else {
				static_assert( ::std::is_member_object_pointer_v<decltype( f )> );
				static_assert( sizeof...( args ) == 0 );
				if constexpr ( ::std::is_base_of_v<T, ::std::decay_t<T1>> )
					return ::std::forward<T1>( t1 ).*f;
				else if constexpr ( is_reference_wrapper_v<::std::decay_t<T1>> )
					return t1.get().*f;
				else
					return ( *::std::forward<T1>( t1 ) ).*f;
			}
		}

		template <class F, class... Args>
		constexpr decltype( auto ) constexpr_invoke_impl( F &&f, Args &&...args ) {
			return ::std::forward<F>( f )( ::std::forward<Args>( args )... );
		}

	} // namespace detail

	/// TODO: Come C++20, ::std::invoke becomes constexpr so retire this
	template <class F, class... Args>
	constexpr ::std::invoke_result_t<F, Args...> constexpr_invoke( F &&f, Args &&...args ) noexcept(
	  ::std::is_nothrow_invocable_v<F, Args...> ) {
		return detail::constexpr_invoke_impl( ::std::forward<F>( f ), ::std::forward<Args>( args )... );
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_CONSTEXPR_INVOKE_HPP
