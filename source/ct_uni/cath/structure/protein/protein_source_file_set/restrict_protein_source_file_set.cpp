/// \file
/// \brief The restrict_protein_source_file_set class definitions

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

#include "restrict_protein_source_file_set.hpp"

#include "cath/chopping/region/region.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::file;

using boost::none;
using std::ostream;
using std::string;

/// \brief For a restrict_protein_source_file_set, implement do_read_files in terms of do_read_and_restrict_files()
protein restrict_protein_source_file_set::do_read_files(const data_file_path_map &prm_filename_of_data_file, ///< The pre-loaded map of file types to filenames
                                                        const string             &prm_protein_name,          ///< The name of the protein that is to be read from files
                                                        ostream                  &prm_stderr                 ///< The ostream to which any warnings/errors should be written
                                                        ) const {
	return do_read_and_restrict_files(
		prm_filename_of_data_file,
		prm_protein_name,
		none,
		prm_stderr
	);
}
