/// \file
/// \brief The type_to_string() header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TYPE_TO_STRING_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TYPE_TO_STRING_HPP

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/core/demangle.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/range/algorithm/mismatch.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/cpp20/make_array.hpp"

#include <array>
#include <string>
#include <typeinfo>

namespace cath {
	namespace common {

		template <typename T>
		std::string type_to_string();

		namespace detail {

			/// \brief A typelist struct to be used in identifying other template wrapper names
			template <typename... Ts>
			struct some_other_class final {};

			/// \brief Tidy up the specified name string to remove common inline implementation namespaces
			inline void tidy_string(std::string &arg_string ///< The string to tidy up
			                        ) {
				if ( boost::algorithm::starts_with( arg_string, "std::__1::" ) ) {
					arg_string = "std::" + arg_string.substr( 10 );
				}
				if ( boost::algorithm::starts_with( arg_string, "std::__debug::" ) ) {
					arg_string = "std::" + arg_string.substr( 14 );
				}
			}

			/// \brief Tidy up a copy of the specified name string to remove common inline implementation namespaces
			inline std::string tidy_string_copy(std::string arg_string ///< The string to tidy up
			                                    ) {
				tidy_string( arg_string );
				return arg_string;
			}

			/// \brief Return a string with the name of the template wrapper T
			///
			/// Note: this requires that T<Ts...> be valid
			template <template <typename...> class T,
			          typename... Ts>
			std::string template_wrapper_name() {
				const std::string soc_template_fullname = ::boost::core::demangle( typeid( some_other_class<Ts...> ).name() );
				const std::string new_template_fullname = ::boost::core::demangle( typeid( T               <Ts...> ).name() );

				const auto mistmatch_rev_itr = boost::range::mismatch(
					soc_template_fullname | boost::adaptors::reversed,
					new_template_fullname | boost::adaptors::reversed
				).second;

				return tidy_string_copy( std::string{
					common::cbegin( new_template_fullname ),
					mistmatch_rev_itr.base()
				} );
			}

			/// \brief Primary template for making the name of the specified type
			template <typename T>
			struct type_to_string_impl final {
				std::string operator()() const {
					return tidy_string_copy( ::boost::core::demangle( typeid( T ).name() ) );
				}
			};

			/// \brief Specialisation for making the name of a type that's a template of types 
			template <template <typename...> class T,
			          typename... Ts>
			struct type_to_string_impl<T<Ts...>> final {
				std::string operator()() const {
					const auto wrapper_name = template_wrapper_name<T, Ts...>();
					const auto param_names = make_array( type_to_string<Ts>()... );
					return wrapper_name
						+ "<"
						+ boost::algorithm::join(
							param_names,
							", "
						)
						+ ">";
				}
			};
		} // namespace detail

		/// \brief Generate a string containing a name of the type T
		template <typename T>
		inline std::string type_to_string() {
			return detail::type_to_string_impl<T>{}();
		}

	} // namespace common
} // namespace cath

#endif
