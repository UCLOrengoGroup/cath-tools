/// \file
/// \brief The tribool functions header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_TRIBOOL_TRIBOOL_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_TRIBOOL_TRIBOOL_HPP

#include <boost/logic/tribool.hpp>

namespace cath {
	namespace common {

		/// \brief Return whether the specified tribool is true
		///
		/// \TODO Come Boost >= 1.58.0 (even on Travis-CI), make this constexpr
		inline bool is_true(const boost::logic::tribool &arg_tribool ///< The tribool to query
		                    ) {
			return static_cast<bool>(   arg_tribool );
		}

		/// \brief Return whether the specified tribool is false
		///
		/// \TODO Come Boost >= 1.58.0 (even on Travis-CI), make this constexpr
		inline bool is_false(const boost::logic::tribool &arg_tribool ///< The tribool to query
		                     ) {
			return static_cast<bool>( ! arg_tribool );
		}

		/// \brief Return whether the specified tribool is not true
		///
		/// \TODO Come Boost >= 1.58.0 (even on Travis-CI), make this constexpr
		inline bool is_not_true(const boost::logic::tribool &arg_tribool ///< The tribool to query
		                        ) {
			return ! is_true( arg_tribool );
		}

		/// \brief Return whether the specified tribool is not false
		///
		/// \TODO Come Boost >= 1.58.0 (even on Travis-CI), make this constexpr
		inline bool is_not_false(const boost::logic::tribool &arg_tribool ///< The tribool to query
		                         ) {
			return ! is_false( arg_tribool );
		}

		/// \brief Return whether the specified tribool is not indeterminate
		///
		/// \TODO Come Boost >= 1.58.0 (even on Travis-CI), make this constexpr
		inline bool is_not_indeterminate(const boost::logic::tribool &arg_tribool ///< The tribool to query
		                                 ) {
			return ! boost::logic::indeterminate( arg_tribool );
		}

	} // namespace common
} // namespace cath

#endif
