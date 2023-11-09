/// \file
/// \brief The make_ptree_of header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_MAKE_PTREE_OF_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_MAKE_PTREE_OF_HPP

#include <boost/property_tree/ptree.hpp>

namespace cath::common {

	/// \brief Make a ptree of the specified value
	///
	/// \tparam T must have an associated `save_to_ptree(ptree &, const T &)` non-member function
	template <typename T>
	boost::property_tree::ptree make_ptree_of(const T &prm_val ///< The value to represent in the ptree
	                                          ) {
		boost::property_tree::ptree new_ptree;
		save_to_ptree( new_ptree, prm_val );
		return new_ptree;
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_MAKE_PTREE_OF_HPP
