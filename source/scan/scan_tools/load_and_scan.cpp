/// \file
/// \brief The load_and_scan class definitions

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

#include "load_and_scan.hpp"

#include "scan/scan_action/record_scores_scan_action.hpp"
#include "scan/scan_tools/load_and_scan_metrics.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::scan;
using namespace std;

/// \brief TODOCUMENT
void load_and_scan::perform_load() {
	stringstream stderr_ostream;
	const auto query_load_result = query_protein_loader.load_proteins( stderr_ostream );
	const auto match_load_result = match_protein_loader.load_proteins( stderr_ostream );
	query_proteins = query_load_result.first;
	match_proteins = match_load_result.first;
	load_files_duration = query_load_result.second + match_load_result.second;
}

/// \brief TODOCUMENT
void load_and_scan::perform_scan() {
	const auto scan_results = scan_ptr->perform_scan( get_query_proteins(), get_match_proteins() );
	the_scan_metrics = scan_results.second;
}

/// \brief TODOCUMENT
load_and_scan::load_and_scan(const protein_list_loader &arg_query_protein_list_loader, ///< TODOCUMENT
                             const protein_list_loader &arg_match_protein_list_loader, ///< TODOCUMENT
                             const scan_type           &arg_scan_type                  ///< TODOCUMENT
                             ) : query_protein_loader( arg_query_protein_list_loader ),
                                 match_protein_loader( arg_match_protein_list_loader ),
                                 scan_ptr            ( arg_scan_type.clone()         ) {
	perform_load();
	perform_scan();
}

/// \brief TODOCUMENT wireframe quads unbind torsional taylor constness
const protein_list & load_and_scan::get_query_proteins() const {
	return *query_proteins;
}

/// \brief TODOCUMENT
const protein_list & load_and_scan::get_match_proteins() const {
	return *match_proteins;
}

/// \brief TODOCUMENT
load_and_scan_metrics load_and_scan::get_load_and_scan_metrics() const {
	return {
		*load_files_duration,
		*the_scan_metrics
	};
}
