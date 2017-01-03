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

#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/protein_source_file_set/protein_source_file_set.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <chrono>

using namespace boost::filesystem;
using namespace cath;
using namespace std;

/// \brief TODOCUMENT
protein_list_loader::protein_list_loader(const protein_source_file_set &arg_source_file_set, ///< TODOCUMENT
                                         const path                    &arg_data_dir,        ///< TODOCUMENT
                                         const str_vec                 &arg_protein_names    ///< TODOCUMENT
                                         ) : source_file_set_ptr ( arg_source_file_set.clone() ),
                                             data_dir            ( arg_data_dir                ),
                                             protein_names       ( arg_protein_names           ) {
}

/// \brief TODOCUMENT
pair<protein_list, hrc_duration> protein_list_loader::load_proteins(ostream &arg_stderr ///< TODOCUMENT
                                                                    ) const {
	const auto scan_starttime = std::chrono::high_resolution_clock::now();
	const auto the_proteins   = read_proteins_from_files(
		*source_file_set_ptr,
		data_dir,
		protein_names,
		ref( arg_stderr )
	);
	return make_pair(
		the_proteins,
		std::chrono::high_resolution_clock::now() - scan_starttime
	);
}

