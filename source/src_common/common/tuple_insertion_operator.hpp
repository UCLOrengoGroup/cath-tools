/// \file
/// \brief The tuple insertion operator header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TUPLE_INSERTION_OPERATOR_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TUPLE_INSERTION_OPERATOR_HPP

#include <boost/algorithm/string/join.hpp>
#include <boost/core/demangle.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <utility>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Implementation function to generate a string with a comma-separated list of the tuple's values
			template <typename Tpl, size_t... Index>
			std::string tuple_values_to_string_impl(const Tpl &arg_tuple,         ///< The tuple to be described
			                                        std::index_sequence<Index...> ///< An index_sequence matching the indices of Tpl
			                                        ) {
				using std::to_string;
				const auto a = { std::to_string( std::get<Index>( arg_tuple ) )... };
				return boost::algorithm::join( a, ", " );
			}

		} // namespace detail


		/// \brief Generate a string describing the specified tuple
		///
		/// Example output: "std::tuple<m, d>(3, 2.300000)"
		///
		/// \todo Come C++17, use fold expressions to simplify this
		template <typename... Ts>
		std::string tuple_to_string(const std::tuple<Ts...> &arg_tuple ///< The tuple to describe
		                            ) {
			const auto type_names = { ::boost::core::demangle( typeid( Ts ).name() )... };

			return
				  "std::tuple<"
				+ boost::algorithm::join(
					type_names,
					", "
				)
				+ ">("
				+ detail::tuple_values_to_string_impl(
					arg_tuple,
					std::make_index_sequence< sizeof...( Ts ) >{}
				)
				+ ")";
		}

	} // namespace common
} // namespace cath

namespace std {

	/// \brief An insertion operator for a tuple
	///
	/// It's rather unpleasant to be putting this into the std namespace, which is bad practice
	/// Unfortunately, this is the way to get things to work in the Boost Test library ATM
	///
	/// \TODO Come Boost 1.64.0, use `boost_test_print_type` instead, see:
	///  * http://www.boost.org/libs/test/doc/html/boost_test/test_output/test_tools_support_for_logging/testing_tool_output_disable.html#boost_test.test_output.test_tools_support_for_logging.testing_tool_output_disable.user_type_customization_point_fo
	///  * https://svn.boost.org/trac10/ticket/12540
	///
	/// \TODO Move this into a module for test code - this shouldn't be included with any non-test code
	template <typename... Ts>
	ostream & operator<<(ostream            &arg_os,
	                     const tuple<Ts...> &arg_tuple
	                     ) {
		arg_os << cath::common::tuple_to_string( arg_tuple  );
		return arg_os;
	}

} // namespace std

#endif
