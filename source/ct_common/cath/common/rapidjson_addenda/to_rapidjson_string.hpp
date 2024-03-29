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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_TO_RAPIDJSON_STRING_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_TO_RAPIDJSON_STRING_HPP

#include "cath/common/rapidjson_addenda/string_of_rapidjson_write.hpp"

namespace cath::common {

	// Declaration of function template required here to help to_rapidjson_string find the partial
	// specialisation (because there are issues with function templates of values where the compiler
	// can't do the lookup until it's known to be a function but can't know that until its done the lookup)
	template <json_style Style, typename T>
	void write_to_rapidjson(rapidjson_writer<Style> &,
	                        const T & = std::declval<std::enable_if<false, T>>() );

	namespace detail {

		/// \brief Implementation of to_rapidjson_string
		template <json_style Style, typename T, typename... As>
		std::string to_rapidjson_string_impl(const T       &   prm_value,       ///< The value to write
		                                     const size_t  &   prm_extra_depth, ///< The number of levels of depth
		                                     As           &&...prm_args         ///< Any further arguments
		                                     ) {
			return string_of_rapidjson_write<Style>(
				[&] (rapidjson_writer<Style> &the_writer) {
					write_to_rapidjson(
						the_writer,
						prm_value,
						std::forward< As >( prm_args )...
					);
				},
				prm_extra_depth
			);
		}

	} // namespace detail

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
	/// This overload is for providing a default prm_extra_depth when no extra arguments are specified
	///
	/// \tparam T must be a type for which a write_to_rapidjson has been partially specialised
	template <json_style Style, typename T>
	std::string to_rapidjson_string(const T       &prm_value,          ///< The value to write
	                                const size_t  &prm_extra_depth = 0 ///< The number of levels of depth
	                                ) {
		return detail::to_rapidjson_string_impl<Style>(
			prm_value,
			prm_extra_depth
		);
	}

	/// \brief Write the specified value to a string using a rapidjson Writer corresponding to json_style Style
	///
	/// This overload is for taking one or more extra arguments to pass to write_to_rapidjson()
	///
	/// \tparam T must be a type for which a write_to_rapidjson has been partially specialised
	template <json_style Style, typename T, typename A, typename... As>
	std::string to_rapidjson_string(const T       &   prm_value,       ///< The value to write
	                                const size_t  &   prm_extra_depth, ///< The number of levels of depth
	                                A            &&   prm_arg,         ///< The first extra argument
	                                As           &&...prm_args         ///< Any further extra arguments
	                                ) {
		return detail::to_rapidjson_string_impl<Style>(
			prm_value,
			prm_extra_depth,
			std::forward< A  >( prm_arg ),
			std::forward< As >( prm_args )...
		);
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_RAPIDJSON_ADDENDA_TO_RAPIDJSON_STRING_HPP
