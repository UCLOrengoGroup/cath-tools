/// \file
/// \brief The lexical_cast_line header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_LEXICAL_CAST_LINE_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_LEXICAL_CAST_LINE_HPP

#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

#include <string>

/// \brief TODOCUMENT
template <typename T>
T lexical_cast_line(std::istream &prm_istream ///< TODOCUMENT
                    ) {
	std::string line;
	getline(prm_istream, line);
	return boost::lexical_cast<T>(line);
}

/// \brief TODOCUMENT
template <typename T>
T lexical_cast_trimmed_line(std::istream &prm_istream ///< TODOCUMENT
                            ) {
	std::string line;
	getline(prm_istream, line);
	return boost::lexical_cast<T>(boost::algorithm::trim_copy(line));
}

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_LEXICAL_CAST_LINE_HPP
