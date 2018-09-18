/// \file
/// \brief The lexical_cast_stringstream class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_LEXICAL_CAST_STRINGSTREAM_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_LEXICAL_CAST_STRINGSTREAM_HPP

#include <sstream>

namespace cath {
	namespace common {

		/// Commenting this out to see if it can be adequately replaced by std::stod()
		///
		/// \todo Check if this has been adequately replaced
		//
		// /// \brief Do a similar job to boost::lexical_cast<> but using a std::stringstream
		// ///
		// /// lexical_cast<double>() seems to (sometimes?) have problems when run under valgrind
		// /// in that it converts "1.000" to a number that's a tiny bit bigger than 1
		// /// (visible using std::setprecision(50))
		// ///
		// /// Might relate to this: https://bugzilla.redhat.com/show_bug.cgi?id=837650 ?
		// template <typename T, typename F>
		// T lexical_cast_stringstream(const F &prm_from ///< TODOCUMENT
		//                             ) {
		// 	std::stringstream the_ss;
		// 	the_ss << prm_from;
		// 	T result;
		// 	the_ss >> result;
		// 	return result;
		// }

	} // namespace common
} // namespace cath

#endif
