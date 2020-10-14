/// \file
/// \brief The protein_from_pdb_and_calc class definitions

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

#include "protein_from_pdb_and_calc.hpp"

#include <boost/filesystem/path.hpp>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/file/data_file.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/protein_io.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;

using ::std::ostream;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method
unique_ptr<protein_source_file_set> protein_from_pdb_and_calc::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Return that this policy requires a PDB file
data_file_vec protein_from_pdb_and_calc::do_get_file_set() const {
	return { data_file::PDB };
}

/// \brief Return that this policy's primary file is the PDB file
data_file protein_from_pdb_and_calc::do_get_primary_file() const {
	return data_file::PDB;
}

/// \brief Return that the equivalent protein_file_combn value for this is PDB_DSSP_SEC
protein_file_combn protein_from_pdb_and_calc::do_get_protein_file_combn() const {
	return protein_file_combn::PDB_AND_CALC;
}

/// \brief Return whether this policy makes proteins that are SSAP-ready (with data loaded for sec, phi/psi accessibility etc)
bool protein_from_pdb_and_calc::do_makes_ssap_ready_protein() const {
	return true;
}

/// \brief Grab the specified PDB filename and then use it in read_protein_from_pdb_and_calc_dssp_and_sec()
protein protein_from_pdb_and_calc::do_read_and_restrict_files(const data_file_path_map &prm_filename_of_data_file, ///< The pre-loaded map of file types to filenames
                                                              const string             &prm_protein_name,          ///< The name of the structure to be loaded
                                                              const region_vec_opt     &prm_regions,               ///< The regions to which the resulting protein should be restricted
                                                              ostream                  &prm_stderr                 ///< The ostream to which warnings/errors should be written
                                                              ) const {
	return read_protein_from_pdb_and_calc_dssp_and_sec(
		prm_filename_of_data_file.at( data_file::PDB ),
		prm_protein_name,
		prm_regions,
		ref( prm_stderr )
	);
}
