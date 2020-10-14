/// \file
/// \brief The scan_type class definitions

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

#include "scan_type.hpp"

#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/scan/scan_action/record_scores_scan_action.hpp"
#include "cath/scan/scan_tools/scan_metrics.hpp"

using namespace cath::common;
using namespace cath::scan;
using namespace std;

/// \brief TODOCUMENT
pair<record_scores_scan_action, scan_metrics> scan_type::perform_scan(const protein_list &prm_query_protein_list, ///< TODOCUMENT
                                                                      const protein_list &prm_match_protein_list  ///< TODOCUMENT
                                                                      ) {
	return do_perform_scan( prm_query_protein_list, prm_match_protein_list );
}

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<scan_type> scan_type::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}
