/// \file
/// \brief The variables_map_contains header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_VARIABLES_MAP_CONTAINS_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_VARIABLES_MAP_CONTAINS_HPP

#include <boost/program_options.hpp>

namespace cath::common {

	/// \brief Return whether the specified key is in the specified variables_map
	inline bool contains(const boost::program_options::variables_map &prm_vm, ///< The variables_map to query
	                     const std::string                           &prm_key ///< The key to search for
	                     ) {
		return ( prm_vm.count( prm_key ) > 0 );
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_PROGRAM_OPTIONS_VARIABLES_MAP_CONTAINS_HPP
