/// \file
/// \brief The back header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_BACK_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_BACK_HPP

#include <boost/range.hpp>
//#include <boost/range/begin.hpp>

namespace cath::common {

	/// \brief Return a non-const reference to the last element of a range
	///
	/// This is a non-member function for ranges that don't provide their own back method.
	/// This is useful for non-member ranges
	///
	/// \todo Come C++14 (in GCC > v4.8), switch to std::rbegin in <iterator>
	template <typename T>
	inline typename boost::range_reference<T>::type back(T &prm_range ///< The range to query
	                                                     ) {
		return * boost::rbegin( prm_range );
	}

	/// \brief Return a const reference to the first element of a range
	///
	/// \copydetails back()
	///
	/// \todo Come C++14 (in GCC > v4.8), switch to std::crbegin in <iterator>
	template <typename T>
	inline typename boost::range_reference<const T>::type back(const T &prm_range ///< The range to query
	                                                           ) {
		return * boost::const_rbegin( prm_range );
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_BACK_HPP
