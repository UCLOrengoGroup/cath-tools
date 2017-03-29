/// \file
/// \brief The to_rapidjson_string header

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

#ifndef _CATH_TOOLS_SOURCE_COMMON_RAPIDJSON_ADDENDA_TO_RAPIDJSON_STRING_H
#define _CATH_TOOLS_SOURCE_COMMON_RAPIDJSON_ADDENDA_TO_RAPIDJSON_STRING_H

#include "common/rapidjson_addenda/string_of_rapidjson_write.hpp"

namespace cath {
	namespace common {

		// Declaration of function template required here to help to_rapidjson_string find the partial
		// specialisation (because there are issues with function templates of values where the compiler
		// can't do the lookup until it's known to be a function but can't know that until its done the lookup)
		template <json_style Style, typename T>
		void write_to_rapidjson(rapidjson_writer<Style> &,
		                        const T & = std::declval<std::enable_if<false, T>>() );

		// // Attempt to compile-time detect whether write_to_rapidjson is (partially) specialised for T
		// // ...was struggling to get this to work (accepts the above primary template for
		// // types for which it hasn't been specialised) and didn't want to waste more time on it.
		// // Not sure it makes sense to test for partial specialisation
		// namespace detail {

		// 	template <typename T,
		// 	          typename = void_t<>
		// 	          >
		// 	struct is_writable_to_rapidjson : std::false_type {};

		// 	template <typename T>
		// 	struct is_writable_to_rapidjson< T, void_t< decltype( write_to_rapidjson<json_style::PRETTY, T>( std::declval<rapidjson_writer<json_style::PRETTY> &>(), std::declval<T>() ) ) > > : std::true_type {};

		// } // namespace detail

		/// \brief Write the specified value to a string using a rapidjson Writer corresponding to json_style Style
		///
		/// \tparam T must be a type for which a write_to_rapidjson has been partially specialised
		template <json_style Style, typename T>
		std::string to_rapidjson_string(const T      &arg_value,          ///< The value to write
		                                const size_t &arg_extra_depth = 0 ///< The number of levels of depth
		                                ) {
			return string_of_rapidjson_write<Style>(
				[&] (rapidjson_writer<Style> &the_writer) { write_to_rapidjson( the_writer, arg_value ); },
				arg_extra_depth
			);
		}


	} // namespace common
} // namespace cath

#endif
