/// \file
/// \brief The res_pair_dirn header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_DIRN_RES_PAIR_DIRN_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_DIRN_RES_PAIR_DIRN_HPP

#include <iosfwd>

namespace cath::scan::detail {

	/// \brief Whether a res_pair's from-residue comes before or after its to-residue
	enum class res_pair_dirn : bool {
		INCREASE, ///< The from-residue comes before the to-residue
		DECREASE  ///< The from-residue comes after  the to-residue
	};

	::std::ostream &operator<<( ::std::ostream &, const res_pair_dirn & );

} // namespace cath::scan::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_DIRN_RES_PAIR_DIRN_HPP
