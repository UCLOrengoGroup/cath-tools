/// \file
/// \brief The read_from_ptree header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_READ_FROM_PTREE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_READ_FROM_PTREE_HPP

#include <boost/property_tree/ptree_fwd.hpp>

namespace cath::common {

	/// \brief TODOCUMENT
	template <typename T>
	T read_from_ptree(const boost::property_tree::ptree &);

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROPERTY_TREE_READ_FROM_PTREE_HPP
