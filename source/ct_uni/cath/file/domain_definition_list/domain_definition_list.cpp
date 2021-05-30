/// \file
/// \brief The domain_definition_list class definitions

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

#include "domain_definition_list.hpp"

#include <filesystem>
#include <fstream>

#include <boost/algorithm/string/trim.hpp>

#include "cath/chopping/chopping_format/sillitoe_chopping_format.hpp"
#include "cath/chopping/domain/domain_definition.hpp"
#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/file/name_set/name_set_list.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/file/pdb/read_domain_def_from_pdb.hpp"

using namespace ::cath;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::boost::algorithm::is_any_of;
using ::boost::algorithm::token_compress_on;
using ::boost::trim;
using ::std::filesystem::path;
using ::std::ifstream;
using ::std::istream;
using ::std::make_pair;
using ::std::nullopt;
using ::std::string;

/// \brief Ctor for domain_definition_list
domain_definition_list::domain_definition_list(domain_definition_vec prm_domain_definitions
                                               ) : domain_definitions { std::move( prm_domain_definitions ) } {
}

/// \brief TODOCUMENT
size_t domain_definition_list::size() const {
	return domain_definitions.size();
}

/// \brief TODOCUMENT
domain_definition_list::const_iterator domain_definition_list::begin() const {
	return cbegin( domain_definitions );
}

/// \brief TODOCUMENT
domain_definition_list::const_iterator domain_definition_list::end() const {
	return cend( domain_definitions );
}

/// \brief TODOCUMENT
domain_definition_list cath::file::parse_domain_definition_file(const path &prm_dom_defn_file ///< TODOCUMENT
                                                                ) {
	ifstream dom_defn_ifstream = open_ifstream( prm_dom_defn_file );
	domain_definition_list new_domain_definition_list = parse_domain_definition_file( dom_defn_ifstream );
	dom_defn_ifstream.close();
	return new_domain_definition_list;
}

/// \brief TODOCUMENT
domain_definition_list cath::file::parse_domain_definition_file(istream &prm_dom_defn_istream ///< TODOCUMENT
                                                                ) {
	string line_string;
	domain_definition_vec domain_definitions;
	while ( getline( prm_dom_defn_istream, line_string ) ) {
		trim( line_string );
		if ( ! line_string.empty() ) {
			const auto line_parts = split_build<str_vec>( line_string, is_any_of( " " ), token_compress_on );
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
pdb_list_name_set_list_pair cath::file::read_domains_from_pdbs(const domain_definition_list &prm_domain_definition_list, ///< TODOCUMENT
                                                               const data_dirs_spec         &prm_data_dirs_spec          ///< TODOCUMENT
                                                               ) {
	const size_t num_domain_definitions = prm_domain_definition_list.size();
	pdb_list pdbs;
	name_set_vec names;
	pdbs.reserve ( num_domain_definitions );
	names.reserve( num_domain_definitions );

	for (const domain_definition &domain_defn : prm_domain_definition_list) {
		const domain &the_domain = domain_defn.get_domain();
		if ( ! has_domain_id( the_domain ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Domain definitions to be read from PDBs do not have domain IDs"));
		}
		auto [file, the_pdb] = read_domain_from_pdb( domain_defn, prm_data_dirs_spec );
		names.emplace_back( std::move( file ), nullopt, get_domain_id( the_domain ) );
		pdbs.push_back( std::move( the_pdb ) );
	}
	return make_pair( pdbs, name_set_list{ std::move( names ) } );
}
