/// \file
/// \brief The log_scan_action class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SCAN_ACTION_LOG_SCAN_ACTION_H
#define _CATH_TOOLS_SOURCE_SCAN_SCAN_ACTION_LOG_SCAN_ACTION_H

#include "scan/detail/res_pair/single_struc_res_pair.h"

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		template <typename T> class TD;
		struct log_scan_action final {
			long long unsigned int num_matches = 0;
			void operator()(const detail::single_struc_res_pair &/*arg_res_pair_a*/,  ///< TODOCUMENT
			                const detail::single_struc_res_pair &/*arg_res_pair_b*/,  ///< TODOCUMENT
			                const index_type                    &/*arg_structure_a*/, ///< TODOCUMENT
			                const index_type                    &/*arg_structure_b*/  ///< TODOCUMENT
			                ) {
				++num_matches;
//				std::cerr << "For query residue pair [query_structure: ";
//				std::cerr << arg_structure_a;
//				std::cerr << "; ";
//				std::cerr << arg_res_pair_a;
//				std::cerr << "], found match [match_structure: ";
//				std::cerr << arg_structure_b;
//				std::cerr << "; ";
//				std::cerr << arg_res_pair_b;
//				std::cerr << "]\n";
			}
		};


	} // namespace scan
} // namespace cath

#endif
