/// \file
/// \brief The pymol_file_superposition_outputter class definitions

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

#include "pymol_file_superposition_outputter.h"

#include "common/clone/make_uptr_clone.h"
#include "common/file/open_fstream.h"
#include "display/viewer/pymol_viewer.h"
#include "display/display_colour/display_colour_list.h"

#include <fstream>

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pymol_file_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pymol_file_superposition_outputter::do_output_superposition(const superposition_context &arg_superposition_context, ///< TODOCUMENT
                                                                 ostream                     &/*arg_ostream*/            ///< TODOCUMENT
                                                                 ) const {
	ofstream pymol_file_ostream;
	open_ofstream(pymol_file_ostream, output_file);

	output_superposition_to_viewer(
		pymol_file_ostream,
		pymol_viewer(),
		the_display_spec,
		arg_superposition_context
	);
	pymol_file_ostream << flush;
	pymol_file_ostream.close();
}

/// \brief TODOCUMENT
bool pymol_file_superposition_outputter::do_involves_display_spec() const {
	return true;
}

/// \brief Ctor for pymol_file_superposition_outputter.
pymol_file_superposition_outputter::pymol_file_superposition_outputter(const path         &arg_output_file, ///< TODOCUMENT
                                                                       const display_spec &arg_display_spec ///< TODOCUMENT
                                                                       ) : output_file      ( arg_output_file  ),
                                                                           the_display_spec ( arg_display_spec ) {
}

