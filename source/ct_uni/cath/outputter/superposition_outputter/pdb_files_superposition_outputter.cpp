/// \file
/// \brief The pdb_files_superposition_outputter class definitions

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

#include "pdb_files_superposition_outputter.hpp"


#include "cath/chopping/region/region.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/superposition/io/superposition_io.hpp"
#include "cath/superposition/superposition_context.hpp"

using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;
using namespace ::cath::sup;

using ::boost::filesystem::path;
using ::boost::string_ref;
using ::std::ostream;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<superposition_outputter> pdb_files_superposition_outputter::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
void pdb_files_superposition_outputter::do_output_superposition(const superposition_context &prm_supn_context, ///< TODOCUMENT
                                                                ostream                     &/*prm_ostream*/   ///< TODOCUMENT
                                                                ) const {
	const pdb_list       pdbs      = get_supn_content_pdbs( prm_supn_context, content_spec );
	const name_set_list &name_sets = get_name_sets( prm_supn_context );
	const str_vec       &names     = get_supn_pdb_file_names( name_sets );

	for (const size_t &pdb_ctr : indices( pdbs.size() ) ) {
		write_superposed_pdb_to_file(
			prm_supn_context.get_superposition(),
			( output_dir / names[ pdb_ctr ] ).string(),
			pdbs[ pdb_ctr ],
			pdb_ctr
		);
	}
}

/// \brief TODOCUMENT
bool pdb_files_superposition_outputter::do_involves_display_spec() const {
	return false;
}

/// \brief Getter for the name of this superposition_outputter
string pdb_files_superposition_outputter::do_get_name() const {
	return "pdb_files_superposition_outputter";
}

/// \brief Ctor for pdb_files_superposition_outputter
pdb_files_superposition_outputter::pdb_files_superposition_outputter(const path                 &prm_output_dir,  ///< TODOCUMENT
                                                                     superposition_content_spec  prm_content_spec ///< The specification of what should be included in the superposition
                                                                     ) : output_dir  { prm_output_dir                },
                                                                         content_spec{ std::move( prm_content_spec ) } {
}

