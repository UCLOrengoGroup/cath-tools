/// \file
/// \brief The pdb_file_superposition_outputter class definitions

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

#include "pdb_file_superposition_outputter.h"

#include "common/clone/make_uptr_clone.h"
#include "common/file/open_fstream.h"
#include "options/outputter/superposition_outputter/ostream_superposition_outputter.h"

#include <fstream>

using namespace boost::filesystem;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pdb_file_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pdb_file_superposition_outputter::do_output_superposition(const superposition_context &arg_superposition_context, ///< TODOCUMENT
                                                               ostream                     &/*arg_ostream*/ ///< TODOCUMENT
                                                               ) const {
	ofstream pdb_ostream;
	open_ofstream(pdb_ostream, output_file);

	ostream_superposition_outputter stream_outputter;
	stream_outputter.output_superposition(arg_superposition_context, pdb_ostream);

	pdb_ostream << flush;
	pdb_ostream.close();
}

/// \brief TODOCUMENT
bool pdb_file_superposition_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Ctor for pdb_file_superposition_outputter.
pdb_file_superposition_outputter::pdb_file_superposition_outputter(const path &arg_output_file
                                                                   ) : output_file(arg_output_file) {
}
