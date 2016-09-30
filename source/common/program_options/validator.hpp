/// \file
/// \brief The validator class header

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

#ifndef VALIDATOR_H_INCLUDED
#define VALIDATOR_H_INCLUDED

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "common/type_aliases.h"

#include <string>
#include <vector>

namespace cath {
	namespace common {

		/// \brief Contain a function for validating lexical_cast-able types
		template <typename T, typename U = T>
		class lex_castable_validator final {
		public:
			static boost::any perform_validate(const boost::any &,
			                                   const str_vec &);
		};

		/// \brief Validate a lexical_cast-able type by attempting to lexical_cast it and throwing an invalid_option_value
		///        if any exception is thrown
		///
		/// This can cut out a lot of the boiler-plate code for types so the validate fn only need contain:
		///
		/// ~~~~~.cpp
		/// arg_value = lex_castable_validator<my_type>::perform_validate( arg_value, arg_value_strings );
		/// ~~~~~
		template <typename T, typename U>
		boost::any lex_castable_validator<T, U>::perform_validate(const boost::any &arg_prev_value,   ///< The previous value (if any)
		                                                          const str_vec    &arg_value_strings ///< The value strings
		                                                          ) {
			// Standard validate boilerplate:
			//  * Make sure no previous assignment to 'a' was made.
			//  * Extract the first string from 'arg_value_strings'.
			//    (If there is more than one string, it's an error, and exception will be thrown.)
			boost::program_options::validators::check_first_occurrence( arg_prev_value );
			const std::string &value_string = boost::program_options::validators::get_single_string( arg_value_strings );

			// Attempt to lexical_cast value_string and if it fails, throw an invalid_option_value exception
			try {
				return U{ boost::lexical_cast<T>( value_string ) };
			}
			catch (...) {
				BOOST_THROW_EXCEPTION( boost::program_options::invalid_option_value( value_string ) );
			}
		}

	} // namespace common
} // namespace cath

#endif
