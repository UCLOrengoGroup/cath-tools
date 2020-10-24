/// \file
/// \brief The invoke header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_INVOKE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_INVOKE_HPP

#include <functional>

namespace cath {
	namespace common {

		// All of this is just copied/pasted from http://en.cppreference.com/w/cpp/utility/functional/invoke
		// (with minor adjustments to indentation)

		// Note that invoke isn't constexpr (though implementers may make their invoke constexpr)
		// This appears to be a complicated issue, eg see conversation below
		// https://www.reddit.com/r/cpp/comments/35g7f6/c17_progress_update/cr44el9/

		namespace detail {
			template <class T>
			struct is_reference_wrapper : std::false_type {};
			template <class U>
			struct is_reference_wrapper<std::reference_wrapper<U>> : std::true_type {};
			// template <class T>
			// constexpr bool is_reference_wrapper_v = is_reference_wrapper<T>::value;

			template <class Base, class T, class Derived, class... Args>
			auto INVOKE(T Base::*pmf, Derived&& ref, Args&&... args)
			    noexcept(noexcept((std::forward<Derived>(ref).*pmf)(std::forward<Args>(args)...)))
			 -> std::enable_if_t<std::is_function<T>::value &&
			                     std::is_base_of<Base, std::decay_t<Derived>>::value,
			    decltype((std::forward<Derived>(ref).*pmf)(std::forward<Args>(args)...))>
			{
			      return (std::forward<Derived>(ref).*pmf)(std::forward<Args>(args)...);
			}

			template <class Base, class T, class RefWrap, class... Args>
			auto INVOKE(T Base::*pmf, RefWrap&& ref, Args&&... args)
			    noexcept(noexcept((ref.get().*pmf)(std::forward<Args>(args)...)))
			 -> std::enable_if_t<std::is_function<T>::value &&
			                     is_reference_wrapper<std::decay_t<RefWrap>>::value,
			    decltype((ref.get().*pmf)(std::forward<Args>(args)...))>
			{
			      return (ref.get().*pmf)(std::forward<Args>(args)...);
			}

			template <class Base, class T, class Pointer, class... Args>
			auto INVOKE(T Base::*pmf, Pointer&& ptr, Args&&... args)
			    noexcept(noexcept(((*std::forward<Pointer>(ptr)).*pmf)(std::forward<Args>(args)...)))
			 -> std::enable_if_t<std::is_function<T>::value &&
			                     !is_reference_wrapper<std::decay_t<Pointer>>::value &&
			                     !std::is_base_of<Base, std::decay_t<Pointer>>::value,
			    decltype(((*std::forward<Pointer>(ptr)).*pmf)(std::forward<Args>(args)...))>
			{
			      return ((*std::forward<Pointer>(ptr)).*pmf)(std::forward<Args>(args)...);
			}

			template <class Base, class T, class Derived>
			auto INVOKE(T Base::*pmd, Derived&& ref)
			    noexcept(noexcept(std::forward<Derived>(ref).*pmd))
			 -> std::enable_if_t<!std::is_function<T>::value &&
			                     std::is_base_of<Base, std::decay_t<Derived>>::value,
			    decltype(std::forward<Derived>(ref).*pmd)>
			{
			      return std::forward<Derived>(ref).*pmd;
			}

			template <class Base, class T, class RefWrap>
			auto INVOKE(T Base::*pmd, RefWrap&& ref)
			    noexcept(noexcept(ref.get().*pmd))
			 -> std::enable_if_t<!std::is_function<T>::value &&
			                     is_reference_wrapper<std::decay_t<RefWrap>>::value,
			    decltype(ref.get().*pmd)>
			{
			      return ref.get().*pmd;
			}

			template <class Base, class T, class Pointer>
			auto INVOKE(T Base::*pmd, Pointer&& ptr)
			    noexcept(noexcept((*std::forward<Pointer>(ptr)).*pmd))
			 -> std::enable_if_t<!std::is_function<T>::value &&
			                     !is_reference_wrapper<std::decay_t<Pointer>>::value &&
			                     !std::is_base_of<Base, std::decay_t<Pointer>>::value,
			    decltype((*std::forward<Pointer>(ptr)).*pmd)>
			{
			      return (*std::forward<Pointer>(ptr)).*pmd;
			}

			template <class F, class... Args>
			auto INVOKE(F&& f, Args&&... args)
			    noexcept(noexcept(std::forward<F>(f)(std::forward<Args>(args)...)))
			 -> std::enable_if_t<!std::is_member_pointer<std::decay_t<F>>::value,
			    decltype(std::forward<F>(f)(std::forward<Args>(args)...))>
			{
			      return std::forward<F>(f)(std::forward<Args>(args)...);
			}
		} // namespace detail

		/// \TODO Come C++17, remove this and replace all calls with calls to ::std::invoke
		template< class F, class... ArgTypes >
		auto invoke(F&& f, ArgTypes&&... args)
		    // exception specification for QoI
		    noexcept(noexcept(detail::INVOKE(std::forward<F>(f), std::forward<ArgTypes>(args)...)))
		 -> decltype(detail::INVOKE(std::forward<F>(f), std::forward<ArgTypes>(args)...))
		{
		    return detail::INVOKE(std::forward<F>(f), std::forward<ArgTypes>(args)...);
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP17_INVOKE_HPP
