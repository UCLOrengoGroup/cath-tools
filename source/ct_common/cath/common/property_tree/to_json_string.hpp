/// \file
/// \brief The to_json_string header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_TO_JSON_STRING_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_TO_JSON_STRING_HPP

#include <boost/property_tree/json_parser.hpp>

#include "cath/common/json_style.hpp"
#include "cath/common/property_tree/make_ptree_of.hpp"

#include <string>
#include <sstream>

namespace cath {
	namespace common {

		/// \brief Create a JSON string of the specified value (via a ptree)
		///
		/// \tparam T must have an associated `save_to_ptree(ptree &, const T &)` non-member function
		template <typename T>
		std::string to_json_string(const T          &prm_val,                            ///< The value to represent in the JSON string
		                           const json_style &prm_json_style = json_style::PRETTY ///< The style in which the JSON should be written
		                           ) {
			std::ostringstream json_ss;
			write_json(
				json_ss,
				make_ptree_of( prm_val ),
				( prm_json_style == json_style::PRETTY )
			);
			return json_ss.str();
		}


	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_TO_JSON_STRING_HPP
