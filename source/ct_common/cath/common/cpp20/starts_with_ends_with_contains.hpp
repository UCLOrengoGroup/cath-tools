/// \file
/// \brief The starts_with / ends_with / contains header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_STARTS_WITH_ENDS_WITH_CONTAINS_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_STARTS_WITH_ENDS_WITH_CONTAINS_HPP

#include <string_view>

namespace cath::common {

	constexpr bool starts_with( ::std::string_view x, ::std::string_view y ) noexcept {
		return x.substr( 0, y.size() ) == y;
	}

	constexpr bool ends_with( ::std::string_view x, ::std::string_view y ) noexcept {
		return x.size() >= y.size() && x.compare( x.size() - y.size(), ::std::string_view::npos, y ) == 0;
	}

}

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_CPP20_STARTS_WITH_ENDS_WITH_CONTAINS_HPP
