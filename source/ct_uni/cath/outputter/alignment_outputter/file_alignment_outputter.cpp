/// \file
/// \brief The file_alignment_outputter class definitions

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

#include "file_alignment_outputter.hpp"

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/outputter/alignment_outputter/fasta_ostream_alignment_outputter.hpp"

#include <fstream>

using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::filesystem::path;

/// \brief A standard do_clone method
unique_ptr<alignment_outputter> file_alignment_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void file_alignment_outputter::do_output_alignment(const alignment_context &prm_alignment_context, ///< TODOCUMENT
                                                   ostream                 &/*prm_ostream*/        ///< TODOCUMENT
                                                   ) const {
	ofstream pdb_ostream;
	open_ofstream( pdb_ostream, output_file );

	ostream_alignment_outputter_ptr->output_alignment( prm_alignment_context, pdb_ostream );

	pdb_ostream << flush;
	pdb_ostream.close();
}

/// \brief TODOCUMENT
bool file_alignment_outputter::do_involves_display_spec() const {
	return ostream_alignment_outputter_ptr->involves_display_spec();
}

/// \brief Ctor for file_alignment_outputter.
file_alignment_outputter::file_alignment_outputter(const path                &prm_output_file,        ///< TODOCUMENT
                                                   const alignment_outputter &prm_alignment_outputter ///< TODOCUMENT
                                                   ) : output_file                    ( prm_output_file                 ),
                                                       ostream_alignment_outputter_ptr( prm_alignment_outputter.clone() ) {
}

/// \brief Get a name for this alignment_outputter
string file_alignment_outputter::do_get_name() const {
	return "File[" + ostream_alignment_outputter_ptr->get_name() + "]";
}
