/// \file
/// \brief The set_opt_str_from_prog_opts_try header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_SET_OPT_STR_FROM_PROG_OPTS_TRY_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_SET_OPT_STR_FROM_PROG_OPTS_TRY_HPP

#include <boost/optional.hpp>

#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace common {

		/// \brief If the specified str_opt isn't already set, try the specified function and if any exceptions are thrown
		///        and update the str_opt with a description of the error (using the specified prefix and suffix)
		///
		/// This is a helpful way to provide a consistent approach to different Boost Program Options calls that
		/// may throw an exception to indicate a parsing error
		template <typename Func>
		void set_opt_str_from_prog_opts_try(str_opt           &prm_error_string,  ///< The optional string to update with a description of any errors that occur
		                                    Func             &&prm_function,      ///< The function to perform
		                                    const std::string &prm_prefix_string, ///< The string with which to prefix any error descriptions
		                                    const std::string &prm_suffix_string  ///< The string with which to suffix any error descriptions
		                                    ) {
			if ( ! prm_error_string ) {
				try {
					std::forward<Func>( prm_function )();
				}
				catch (std::exception& e) {
					prm_error_string =
						  prm_prefix_string
						+ e.what()
						+ prm_suffix_string;
				}
				catch (...) {
					prm_error_string =
						  prm_prefix_string
						+ "Caught an unrecognised exception whilst parsing program options."
						+ prm_suffix_string;
				}
			}
		}

	} // namespace common
} // namespace cath

#endif
