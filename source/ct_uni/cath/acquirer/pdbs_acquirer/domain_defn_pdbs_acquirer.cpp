/// \file
/// \brief The domain_defn_pdbs_acquirer class definitions

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

#include "domain_defn_pdbs_acquirer.hpp"

#include <filesystem>
#include <utility>

#include <spdlog/spdlog.h>

#include "cath/chopping/domain/domain_definition.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/file/domain_definition_list/domain_definition_list.hpp"
#include "cath/file/name_set/name_set_list.hpp"
#include "cath/file/options/data_dirs_spec.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"

using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;
using namespace ::std;

using ::std::filesystem::path;

/// \brief A standard do_clone method.
unique_ptr<pdbs_acquirer> domain_defn_pdbs_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief TODOCUMENT
pdb_list_name_set_list_pair domain_defn_pdbs_acquirer::do_get_pdbs_and_names(istream &/*prm_istream*/ ///< TODOCUMENT
                                                                             ) const {
	const domain_definition_list the_dom_defns = parse_domain_definition_file( domain_defn_file );
	::spdlog::warn( "Currently using a hard-coded domain PDB directory : /cath/data/current/pdb" );
	return read_domains_from_pdbs( the_dom_defns, build_data_dirs_spec_of_dir( "/cath/data/current/pdb" ) );
}

/// \brief Ctor for domain_defn_pdbs_acquirer
domain_defn_pdbs_acquirer::domain_defn_pdbs_acquirer(path prm_domain_defn_file ///< TODOCUMENT
                                                     ) : domain_defn_file( std::move( prm_domain_defn_file ) ) {
}

