/// \file
/// \brief The scan_stride class definitions

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

#include "scan_stride.h"

//using namespace cath::scan;
using namespace cath::scan::detail;
//using namespace std;


/// \brief TODOCUMENT
///
/// \relates scan_stride
rep_rep_pair_opt cath::scan::detail::get_from_rep_of_indices(const scan_stride &arg_scan_stride,      ///< TODOCUMENT
                                                             const index_type  &arg_query_from_index, ///< TODOCUMENT
                                                             const index_type  &arg_index_from_index  ///< TODOCUMENT
                                                             ) {
	return get_rep_of_indices(
		arg_scan_stride.get_query_from_strider(),
		arg_query_from_index,
		arg_scan_stride.get_index_from_strider(),
		arg_index_from_index
	);
}

/// \brief TODOCUMENT
///
/// \relates scan_stride
rep_rep_pair_opt cath::scan::detail::get_to_rep_of_indices(const scan_stride &arg_scan_stride,    ///< TODOCUMENT
                                                           const index_type  &arg_query_to_index, ///< TODOCUMENT
                                                           const index_type  &arg_index_to_index  ///< TODOCUMENT
                                                           ) {
	return get_rep_of_indices(
		arg_scan_stride.get_query_to_strider(),
		arg_query_to_index,
		arg_scan_stride.get_index_to_strider(),
		arg_index_to_index
	);
}

