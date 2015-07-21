/// \file
/// \brief The domain_definition_list class definitions

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

#include "domain_definition_list.h"

#include <boost/algorithm/string/trim.hpp>

#include "chopping/chopping_format/sillitoe_chopping_format.h"
#include "chopping/domain/domain_definition.h"
#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/c++14/cbegin_cend.h"
#include "common/file/open_fstream.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_list.h"
#include "file/pdb/pdb_residue.h"

#include <fstream>

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::file;
using namespace cath::opts;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::token_compress_on;

/// \brief Ctor for domain_definition_list
domain_definition_list::domain_definition_list(const domain_definition_vec &arg_domain_definitions
                                               ) : domain_definitions( arg_domain_definitions ) {
}

/// \brief TODOCUMENT
size_t domain_definition_list::size() const {
	return domain_definitions.size();
}

/// \brief TODOCUMENT
domain_definition_list::const_iterator domain_definition_list::begin() const {
	return common::cbegin( domain_definitions );
}

/// \brief TODOCUMENT
domain_definition_list::const_iterator domain_definition_list::end() const {
	return common::cend( domain_definitions );
}

/// \brief TODOCUMENT
domain_definition_list cath::file::parse_domain_definition_file(const path &arg_dom_defn_file ///< TODOCUMENT
                                                                ) {
	ifstream dom_defn_ifstream;
	open_ifstream( dom_defn_ifstream, arg_dom_defn_file );
	const domain_definition_list new_domain_definition_list = parse_domain_definition_file( dom_defn_ifstream );
	dom_defn_ifstream.close();
	return new_domain_definition_list;
}

/// \brief TODOCUMENT
domain_definition_list cath::file::parse_domain_definition_file(istream &arg_dom_defn_istream ///< TODOCUMENT
                                                                ) {
	string line_string;
	domain_definition_vec domain_definitions;
	while ( getline( arg_dom_defn_istream, line_string ) ) {
		trim( line_string );
		if ( ! line_string.empty() ) {
			const str_vec line_parts = split_build<str_vec>( line_string, is_any_of( " " ), token_compress_on );
			if ( line_parts.size() != 3 ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse domain definition from non-empty line that doesn't have three parts"));
			}
			const string &domain_id = line_parts[ 0 ];
			const string &pdb_name  = line_parts[ 1 ];
			const string &chopping  = line_parts[ 2 ];
			domain_definitions.push_back( domain_definition(
				parse_domain( sillitoe_chopping_format(), chopping, domain_id ),
				pdb_name
			) );
		}
	}
	return domain_definition_list( domain_definitions );
}

/// \brief TODOCUMENT
pdb_list_str_vec_pair cath::file::read_domains_from_pdbs(const domain_definition_list  &arg_domain_definition_list, ///< TODOCUMENT
                                                         const data_dirs_options_block &arg_data_dirs_options_block ///< TODOCUMENT
                                                         ) {
	const size_t num_domain_definitions = arg_domain_definition_list.size();
	pdb_list pdbs;
	str_vec  names;
	pdbs.reserve ( num_domain_definitions );
	names.reserve( num_domain_definitions );

	for (const domain_definition &domain_defn : arg_domain_definition_list) {
		const domain &the_domain = domain_defn.get_domain();
		if ( ! has_domain_id( the_domain ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Domain definitions to be read from PDBs do not have domain IDs"));
		}
		names.push_back( get_domain_id( the_domain ) );
		pdbs.push_back ( read_domain_from_pdb( domain_defn, arg_data_dirs_options_block ) );
	}
	return make_pair( pdbs, names );
}

