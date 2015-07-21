/// \file
/// \brief The protein_source_from_wolf_and_sec class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "protein_source_from_wolf_and_sec.h"

#include <boost/filesystem/path.hpp>

#include "common/clone/make_uptr_clone.h"
#include "structure/protein/protein.h"
#include "structure/protein/protein_io.h"
#include "structure/protein/protein_source_file_set/protein_file_combn.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

using namespace boost::filesystem;
using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<protein_source_file_set> protein_source_from_wolf_and_sec::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Return that this policy requires a wolf file and a sec file
data_file_vec protein_source_from_wolf_and_sec::do_get_file_set() const {
	return { data_file::WOLF, data_file::SEC };
}

/// \brief Return that the equivalent protein_file_combn value for this is WOLF_SEC
protein_file_combn protein_source_from_wolf_and_sec::do_get_protein_file_combn() const {
	return protein_file_combn::WOLF_SEC;
}

/// \brief Grab the specified WOLF and SEC filenames and then use them in read_wolf_and_sec_files()
protein protein_source_from_wolf_and_sec::do_read_files(const data_file_path_map &arg_filename_of_data_file, ///< The pre-loaded map of file types to filenames
                                                        const string             &arg_protein_name,          ///< The name of the structure to be loaded
                                                        ostream                  &arg_stderr                 ///< The ostream to which warnings/errors should be written
                                                        ) const {
	const path &wolf_file = arg_filename_of_data_file.at( data_file::WOLF );
	const path &sec_file  = arg_filename_of_data_file.at( data_file::SEC  );
	return read_protein_from_wolf_and_sec_files(
		wolf_file,
		sec_file,
		arg_protein_name,
		arg_stderr
	);
}

