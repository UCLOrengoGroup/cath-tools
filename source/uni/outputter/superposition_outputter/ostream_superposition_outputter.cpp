/// \file
/// \brief The ostream_superposition_outputter class definitions

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

#include "ostream_superposition_outputter.hpp"

#include "chopping/region/region.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "superposition/io/superposition_io.hpp"
#include "superposition/superposition_context.hpp"

using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;

using boost::string_ref;
using std::ostream;
using std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> ostream_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void ostream_superposition_outputter::do_output_superposition(const superposition_context &arg_supn_context, ///< TODOCUMENT
                                                              ostream                     &arg_ostream,      ///< TODOCUMENT
                                                              const string_ref            &/*arg_name*/      ///< A name for the superposition (so users of the superposition know what it represents)
                                                              ) const {
	write_superposed_pdbs_to_ostream(
		arg_ostream,
		arg_supn_context.get_superposition(),
		get_supn_content_pdbs( arg_supn_context, content_spec ),
		sup_pdbs_script_policy::LEAVE_RAW_PDBS,
		chain_relabel_policy::RELABEL
	);
}

/// \brief TODOCUMENT
bool ostream_superposition_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Ctor for ostream_superposition_outputter
ostream_superposition_outputter::ostream_superposition_outputter(superposition_content_spec arg_content_spec ///< The specification of what should be included in the superposition
                                                                 ) : content_spec { std::move( arg_content_spec ) } {
}