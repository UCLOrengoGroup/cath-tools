/// \file
/// \brief The clique header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_CLIQUE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_CLIQUE_HPP

#include <cstddef>
#include <vector>

namespace cath {

	/// \brief Clique secondary structure equivalencies
	struct sec_struc_equivalency final {
		/// \brief TODOCUMENT
		size_t prota_ssnum;

		/// \brief Number of secondary structure
		size_t protb_ssnum;

		/// \brief TODOCUMENT
		char prota_start[ 10 ];

		/// \brief TODOCUMENT
		char prota_end[ 10 ];

		/// \brief TODOCUMENT
		char protb_start[ 10 ];

		/// \brief TODOCUMENT
		char protb_end[ 10 ];
	};

	/// \brief Clique file structure
	struct clique final {
		/// \brief TODOCUMENT
		size_t cliquesize;

		/// \brief Start and ends of secondary structures in proteins a and b
		std::vector<sec_struc_equivalency> equivs;
	};

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_CLIQUE_HPP
