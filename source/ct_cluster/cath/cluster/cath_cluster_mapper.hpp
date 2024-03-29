/// \file
/// \brief The cath_cluster_mapper class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CATH_CLUSTER_MAPPER_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CATH_CLUSTER_MAPPER_HPP

#include "cath/common/type_aliases.hpp"
#include "cath/options/executable/parse_sources.hpp"

#include <iostream>

// clang-format off
namespace cath::clust { class clust_mapping_spec; }
namespace cath::clust { class clustmap_input_spec; }
namespace cath::clust { class clustmap_options; }
namespace cath::clust { class clustmap_output_spec; }
// clang-format on

namespace cath::clust {

	void perform_map_clusters(const str_vec &,
	                          std::istream & = std::cin,
	                          std::ostream & = std::cout,
	                          std::ostream & = std::cerr,
	                          const opts::parse_sources & = opts::parse_sources::CMND_ENV_AND_FILE);

	void perform_map_clusters(const clustmap_options &,
	                          std::istream & = std::cin,
	                          std::ostream & = std::cout,
	                          std::ostream & = std::cerr);

	void perform_map_clusters(const clustmap_input_spec &,
	                          const clust_mapping_spec &,
	                          const clustmap_output_spec &,
	                          std::istream & = std::cin,
	                          std::ostream & = std::cout,
	                          std::ostream & = std::cerr);

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CATH_CLUSTER_MAPPER_HPP
