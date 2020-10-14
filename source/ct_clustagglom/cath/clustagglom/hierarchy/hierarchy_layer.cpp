/// \file
/// \brief The hierarchy class definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
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

#include "hierarchy_layer.hpp"

#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/size_t_literal.hpp"

using namespace ::cath;
using namespace ::cath::common;

using ::boost::algorithm::join;
using ::boost::format;
using ::std::ostream;
using ::std::string;

/// \brief Generate string describing each of the groups in the specified hierarchy_layer
///
/// \relates hierarchy_layer
str_vec cath::clust::to_strings(const hierarchy_layer &prm_hierarchy_layer ///< The hierarchy_layer to describe
                                ) {
	using ::std::to_string;
	return transform_build<str_vec>(
		indices( prm_hierarchy_layer.size() ),
		[&] (const size_t &x) {
			return
				  "GROUP "
				+ ( format( R"(%3d)" ) % x ).str()
				+ ": "
				+ to_string( prm_hierarchy_layer[ x ] );
		}
	);
}

/// \brief Generate a string describing the specified hierarchy_layer
///
/// \relates hierarchy_layer
string cath::clust::to_string(const hierarchy_layer &prm_hierarchy_layer, ///< The hierarchy_layer to describe
                              const string          &prm_join_string      ///< The string with which to join the strings describing the groups
                              ) {
	return join(
		to_strings( prm_hierarchy_layer ),
		prm_join_string
	);
}

/// \brief Insert a description of the specified hierarchy_layer into the specified ostream
///
/// \relates hierarchy_layer
ostream & cath::clust::operator<<(ostream               &prm_os,             ///< The ostream into which the description should be inserted
                                  const hierarchy_layer &prm_hierarchy_layer ///< The hierarchy_layer to describe
                                  ) {
	prm_os << to_string( prm_hierarchy_layer );
	return prm_os;
}
