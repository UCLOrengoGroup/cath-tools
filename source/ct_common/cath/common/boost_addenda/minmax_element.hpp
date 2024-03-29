/// \file
/// \brief The minmax_element header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_MINMAX_ELEMENT_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_MINMAX_ELEMENT_HPP

#include <iterator>

#include <boost/algorithm/minmax_element.hpp>

namespace cath::common {

	/// \brief TODOCUMENT
	///
	/// \param prm_range TODOCUMENT
	template <typename RNG, typename FN>
	inline auto minmax_element( RNG &prm_range ) {
		return boost::minmax_element( ::std::cbegin( prm_range ), ::std::cend( prm_range ) );
	}

	/// \brief TODOCUMENT
	///
	/// \param prm_range    TODOCUMENT
	/// \param prm_function TODOCUMENT
	template <typename RNG, typename FN>
	inline auto minmax_element( RNG &prm_range, FN prm_function ) {
		return boost::minmax_element( ::std::cbegin( prm_range ), ::std::cend( prm_range ), prm_function );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_MINMAX_ELEMENT_HPP
