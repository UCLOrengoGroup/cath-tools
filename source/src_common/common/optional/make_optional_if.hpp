/// \file
/// \brief The make_optional_if() header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_OPTIONAL_MAKE_OPTIONAL_IF_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_OPTIONAL_MAKE_OPTIONAL_IF_HPP

#include <boost/optional.hpp>

#include "common/cpp17/invoke.hpp"

namespace cath {
	namespace common {

		/// \brief Make an optional value from a bool and a value
		///
		/// This should be used over boost::make_optional to make it easy to switch
		/// to std::optional (which doesn't provide this bool version of make_optional)
		///
		/// \todo Apply this across the code base
		///
		/// This could be enhanced:
		///  * Allow 0 or more values after the bool and perfect-forward them to the relevant ctor
		///  * Allow the resulting optional's value_type to be specified but have it default to the
		///    (decay_t<> of the) type of the (first) value after the bool
		///  * If combining both of the above, generate a sensible compile time error if 0
		///    arguments are passed and no template parameter is specified
		template <class T>
		inline boost::optional<T> make_optional_if(const bool &arg_condition, ///< TODOCUMENT
		                                           const T    &arg_value      ///< TODOCUMENT
		                                           ) {
			return arg_condition
				? boost::optional<T>{ arg_value }
				: boost::optional<T>{           };
		}

		/// \brief Make an optional from a bool and an nullary invokable, which is only invoked
		///        if the bool value is true
		///
		/// \todo Apply this across the code base
		///
		/// This could be enhanced:
		///  * Allow the resulting optional's value_type to be specified but have it default to the type
		///    of the std::decay_t< std::result_of_t<> > of the Fn
		template <class Fn>
		inline auto make_optional_if_fn(const bool  &arg_condition, ///< TODOCUMENT
		                                Fn         &&arg_fn         ///< TODOCUMENT
		                                ) {
			using return_type = boost::optional< std::decay_t< std::result_of_t< Fn && () > > >;
			return arg_condition
				? return_type{ invoke( std::forward<Fn>( arg_fn ) ) }
				: return_type{                                      };
		}

	} // namespace common
} // namespace cath

#endif
