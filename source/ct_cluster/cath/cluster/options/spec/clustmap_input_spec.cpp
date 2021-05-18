/// \file
/// \brief The clustmap_input_spec class definitions

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

#include "clustmap_input_spec.hpp"

#include "cath/cluster/options/options_block/clustmap_input_options_block.hpp"

using namespace ::cath;
using namespace ::cath::clust;

using ::std::nullopt;

/// \brief Return a string explaining why the specified clustmap_input_spec is invalid or nullopt if it isn't
///
/// \relates clustmap_input_spec
str_opt cath::clust::get_invalid_description(const clustmap_input_spec &prm_clustmap_input_spec ///< The clustmap_input_spec to query
                                             ) {
	if ( prm_clustmap_input_spec.get_map_from_clustmemb_file() && prm_clustmap_input_spec.get_read_batches_from_input() ) {
		return "Cannot specify a map-from cluster-membership file (--"
			+ clustmap_input_options_block::PO_MAP_FROM_CLUSTMEMB_FILE
			+ ") when reading batches from input (--"
			+ clustmap_input_options_block::PO_READ_BATCHES_FROM_INPUT
			+ ")";
	}

	return nullopt;
}