/// \file
/// \brief The from_json_string header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_FROM_JSON_STRING_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_FROM_JSON_STRING_HPP

#include <sstream>
#include <string>

#include <boost/property_tree/json_parser.hpp>

#include "cath/common/property_tree/read_from_ptree.hpp"

namespace cath {
	namespace common {

		/// \brief Build a T from a JSON string (via a ptree)
		///
		/// Requires that there is specialisation of read_from_ptree<> for T
		template <typename T>
		T from_json_string(const std::string &prm_json_string ///< The JSON string from which the T should be read
		                   ) {
			boost::property_tree::ptree tree;
			std::istringstream in_ss( prm_json_string );
			boost::property_tree::json_parser::read_json( in_ss, tree );
			return ::cath::common::read_from_ptree<T>( tree );
		}

	} // namespace common
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_FROM_JSON_STRING_HPP
