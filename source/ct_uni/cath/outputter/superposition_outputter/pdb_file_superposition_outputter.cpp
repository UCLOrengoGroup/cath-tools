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

#include "pdb_file_superposition_outputter.hpp"

#include <filesystem>
#include <fstream>
#include <utility>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/outputter/superposition_outputter/ostream_superposition_outputter.hpp"

using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::boost::string_ref;
using ::std::filesystem::path;
using ::std::flush;
using ::std::ofstream;
using ::std::ostream;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pdb_file_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pdb_file_superposition_outputter::do_output_superposition(const superposition_context &prm_superposition_context, ///< TODOCUMENT
                                                               ostream                     &/*prm_ostream*/            ///< TODOCUMENT
                                                               ) const {
	ofstream pdb_ostream = open_ofstream( output_file );

	ostream_superposition_outputter stream_outputter{ content_spec };
	stream_outputter.output_superposition( prm_superposition_context, pdb_ostream );

	pdb_ostream << flush;
	pdb_ostream.close();
}

/// \brief TODOCUMENT
bool pdb_file_superposition_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Getter for the name of this superposition_outputter
string pdb_file_superposition_outputter::do_get_name() const {
	return "pdb_file_superposition_outputter";
}

/// \brief Ctor for pdb_file_superposition_outputter
pdb_file_superposition_outputter::pdb_file_superposition_outputter(path prm_output_file, ///< TODOCUMENT
                                                                   superposition_content_spec  prm_content_spec ///< The specification of what should be included in the superposition
                                                                   ) : output_file  { std::move( prm_output_file  ) },
                                                                       content_spec { std::move( prm_content_spec ) } {
}
