/// \file
/// \brief The regex_count header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_REGEX_REGEX_COUNT_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_REGEX_REGEX_COUNT_HPP

#include <regex>

namespace cath::common {

	/// \brief Return the number of matches to the specified regex in the specified string
	///
	/// \param prm_string The string to search
	/// \param prm_regex  The regex to count matches for
	inline size_t regex_count( const std::string &prm_string, const std::regex &prm_regex ) {
		return static_cast<size_t>(
		  std::distance( std::sregex_iterator{ ::std::cbegin( prm_string ), ::std::cend( prm_string ), prm_regex },
		                 std::sregex_iterator{} ) );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_REGEX_REGEX_COUNT_HPP
