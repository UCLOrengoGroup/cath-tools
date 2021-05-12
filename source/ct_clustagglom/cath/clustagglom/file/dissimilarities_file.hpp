/// \file
/// \brief The dissimilarities_file class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE_DISSIMILARITIES_FILE_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE_DISSIMILARITIES_FILE_HPP

#include <filesystem>

#include <boost/optional.hpp>

#include "cath/clustagglom/link_dirn.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath { namespace clust { class links; } }
namespace cath { namespace common { class id_of_str_bidirnl; } }

namespace cath {
	namespace clust {

		links parse_dissimilarities(std::istream &,
		                            common::id_of_str_bidirnl &,
		                            const link_dirn &,
		                            const size_t & = 2);

		links parse_dissimilarities(const std::string &,
		                            common::id_of_str_bidirnl &,
		                            const link_dirn &,
		                            const size_t & = 2);

		links parse_dissimilarities(const ::std::filesystem::path &,
		                            common::id_of_str_bidirnl &,
		                            const link_dirn &,
		                            const size_t & = 2);

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_FILE_DISSIMILARITIES_FILE_HPP
