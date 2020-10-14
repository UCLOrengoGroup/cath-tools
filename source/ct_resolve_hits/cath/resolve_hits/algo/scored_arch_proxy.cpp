/// \file
/// \brief The scored_arch_proxy class definitions

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

#include "scored_arch_proxy.hpp"

#include <boost/algorithm/string/join.hpp>

#include "cath/common/boost_addenda/range/adaptor/lexical_casted.hpp"

#include <string>

using namespace ::cath::common;

using ::boost::algorithm::join;
using ::std::ostream;
using ::std::string;

/// \brief Generate a string describing the specified scored_arch_proxy
///
/// \relates scored_arch_proxy
string cath::rslv::to_string(const scored_arch_proxy &prm_scored_arch_proxy ///< The scored_arch_proxy to describe
                             ) {
	return "scored_arch_proxy[hit indices: "
		+ join( prm_scored_arch_proxy | lexical_casted<string>(), ", " )
		+ "]";
}

/// \brief Insert a description of the specified scored_arch_proxy into the specified ostream
///
/// \relates scored_arch_proxy
ostream & cath::rslv::operator<<(ostream                 &prm_os,               ///< The ostream into which the description should be inserted
                                 const scored_arch_proxy &prm_scored_arch_proxy ///< The scored_arch_proxy to describe
                                 ) {
	prm_os << to_string( prm_scored_arch_proxy );
	return prm_os;
}
