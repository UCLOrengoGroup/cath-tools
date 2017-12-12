/// \file
/// \brief The make_string_ref header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_MAKE_STRING_REF_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_MAKE_STRING_REF_H

#include <boost/utility/string_ref.hpp>

#include "common/debug_numeric_cast.hpp"

namespace cath {
	namespace common {

		/// \brief Make a string_ref from the specified begin and end string iterators
		inline boost::string_ref make_string_ref(const std::string::const_iterator &arg_begin, ///< A string iterator the start              of the region of string to which the string_ref should refer
		                                         const std::string::const_iterator &arg_end    ///< A string iterator the end (one-past-end) of the region of string to which the string_ref should refer
		                                         ) {
			return {
				&*arg_begin,
				debug_numeric_cast<size_t>( std::distance( arg_begin, arg_end ) )
			};
		}

		/// \brief Make a string_ref from the specified pair of begin and end string iterators
		inline boost::string_ref make_string_ref(const std::pair<std::string::const_iterator, std::string::const_iterator> &arg_itrs ///< String iterators to the start and end (one-past-end) the region of string to which the string_ref should refer
		                                         ) {
			return make_string_ref( arg_itrs.first, arg_itrs.second );
		}

		/// \brief Make a string_ref from the specified begin and end string_ref iterators
		inline boost::string_ref make_string_ref(const boost::string_ref::const_iterator &arg_begin, ///< A string_ref iterator the start              of the region of string_ref to which the string_ref should refer
		                                         const boost::string_ref::const_iterator &arg_end    ///< A string_ref iterator the end (one-past-end) of the region of string_ref to which the string_ref should refer
		                                         ) {
			return {
				&*arg_begin,
				debug_numeric_cast<size_t>( std::distance( arg_begin, arg_end ) )
			};
		}

	} // namespace common
} // namespace cath

#endif
