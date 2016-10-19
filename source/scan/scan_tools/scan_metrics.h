/// \file
/// \brief The scan_metrics class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SCAN_TOOLS_SCAN_METRICS_H
#define _CATH_TOOLS_SOURCE_SCAN_SCAN_TOOLS_SCAN_METRICS_H

#include "common/chrono/chrono_type_aliases.h"
#include "scan/detail/scan_type_aliases.h"

#include <map>

namespace cath {
	namespace scan {

			/// \brief TODOCUMENT
			class scan_metrics final {
			private:
				/// \brief TODOCUMENT
				enum class scan_build_type {
					QUERY_STRUCS, ///< TODOCUMENT
					QUERY_INDEX,  ///< TODOCUMENT
					INDEX_STRUCS, ///< TODOCUMENT
					INDEX_INDEX   ///< TODOCUMENT
				};

				/// \brief TODOCUMENT
				std::map<scan_build_type, durn_mem_pair> build_durns_and_sizes;

				/// \brief TODOCUMENT
				hrc_duration scan_durn;

				/// \brief TODOCUMENT
				// hrc_duration_opt  align_all_durn;

				const durn_mem_pair & get_build_durn_and_size(const scan_build_type &) const;

			public:
				scan_metrics(const durn_mem_pair &,
				             const durn_mem_pair &,
				             const durn_mem_pair &,
				             const durn_mem_pair &,
				             const hrc_duration &);

				const durn_mem_pair & get_query_strucs_metrics() const;
				const durn_mem_pair & get_query_index_metrics() const;
				const durn_mem_pair & get_index_strucs_metrics() const;
				const durn_mem_pair & get_index_index_metrics() const;
				const hrc_duration & get_scan_durn() const;
			};

	} // namespace scan
} // namespace cath

#endif
