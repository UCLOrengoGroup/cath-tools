/// \file
/// \brief The protein_list_loader class definitions

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

#include "protein_list_loader.hpp"

#include <chrono>
#include <filesystem>

#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_list.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_source_file_set.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::std;

using ::std::filesystem::path;

/// \brief TODOCUMENT
protein_list_loader::protein_list_loader(const protein_source_file_set &prm_source_file_set, ///< TODOCUMENT
                                         const path                    &prm_data_dir,        ///< TODOCUMENT
                                         str_vec                        prm_protein_names    ///< TODOCUMENT
                                         ) : source_file_set_ptr ( prm_source_file_set.clone()    ),
                                             data_dir            ( prm_data_dir                   ),
                                             protein_names       ( std::move( prm_protein_names ) ) {
}

/// \brief TODOCUMENT
pair<protein_list, hrc_duration> protein_list_loader::load_proteins(ostream &prm_stderr ///< TODOCUMENT
                                                                    ) const {
	const auto scan_starttime = std::chrono::high_resolution_clock::now();
	const auto the_proteins   = read_proteins_from_files(
		*source_file_set_ptr,
		data_dir,
		protein_names,
		ref( prm_stderr )
	);
	return make_pair(
		the_proteins,
		std::chrono::high_resolution_clock::now() - scan_starttime
	);
}

