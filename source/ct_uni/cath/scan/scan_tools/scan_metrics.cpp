/// \file
/// \brief The scan_metrics class definitions

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

#include "scan_metrics.hpp"

#include <boost/units/quantity.hpp>

// #include "cath/scan/detail/scan_type_aliases.hpp"
// #include "cath/scan/scan_action/record_scores_scan_action.hpp"
// #include "cath/scan/scan_index.hpp"
// #include "cath/scan/scan_policy.hpp"
// #include "cath/scan/scan_query_set.hpp"
// #include "cath/scan/scan_stride.hpp"
// #include "cath/structure/geometry/angle.hpp"

// #include <chrono>

using namespace ::cath;
//using namespace ::cath::common;
using namespace ::cath::scan;
//using namespace ::std;

/// \brief TODOCUMENT
const durn_mem_pair & scan_metrics::get_build_durn_and_size(const scan_build_type &prm_scan_build_type ///< TODOCUMENT
                                                            ) const {
	return build_durns_and_sizes.at( prm_scan_build_type );
}

/// \brief TODOCUMENT
scan_metrics::scan_metrics(const durn_mem_pair &prm_build_query_strucs_metrics, ///< TODOCUMENT
                           const durn_mem_pair &prm_build_query_store_metrics,  ///< TODOCUMENT
                           const durn_mem_pair &prm_build_index_strucs_metrics, ///< TODOCUMENT
                           const durn_mem_pair &prm_build_index_store_metrics,  ///< TODOCUMENT
                           const hrc_duration  &prm_scan_durn                   ///< TODOCUMENT
                           ) : build_durns_and_sizes{ {
                               	{ scan_build_type::QUERY_STRUCS, prm_build_query_strucs_metrics },
                               	{ scan_build_type::QUERY_INDEX,  prm_build_query_store_metrics  },
                               	{ scan_build_type::INDEX_STRUCS, prm_build_index_strucs_metrics },
                               	{ scan_build_type::INDEX_INDEX,  prm_build_index_store_metrics  },
                               } },
                               scan_durn      ( prm_scan_durn       ) {
}

/// \brief TODOCUMENT
const durn_mem_pair & scan_metrics::get_query_strucs_metrics() const {
	return get_build_durn_and_size( scan_build_type::QUERY_STRUCS );
}

/// \brief TODOCUMENT
const durn_mem_pair & scan_metrics::get_query_index_metrics() const {
	return get_build_durn_and_size( scan_build_type::QUERY_INDEX );
}

/// \brief TODOCUMENT
const durn_mem_pair & scan_metrics::get_index_strucs_metrics() const {
	return get_build_durn_and_size( scan_build_type::INDEX_STRUCS );
}

/// \brief TODOCUMENT
const durn_mem_pair & scan_metrics::get_index_index_metrics() const {
	return get_build_durn_and_size( scan_build_type::INDEX_INDEX );
}

/// \brief TODOCUMENT
const hrc_duration & scan_metrics::get_scan_durn() const {
	return scan_durn;
}
