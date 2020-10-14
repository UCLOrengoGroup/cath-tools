/// \file
/// \brief The names_file class definitions

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

#include "names_file.hpp"

#include "cath/common/container/id_of_str_bidirnl.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/string/string_parse_tools.hpp"

#include <fstream>
#include <sstream>
#include <string>

using namespace cath;
using namespace cath::clust;
using namespace cath::common;

using boost::filesystem::path;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::pair;
using std::string;

/// \brief Parse the cluster item names from the specified istream, populating the specified name_ider in the process
doub_vec cath::clust::parse_names(istream           &prm_input,    ///< The istream from which to parse the names
                                  id_of_str_bidirnl &prm_name_ider ///< The name_ider to populate from the names data
                                  ) {
	doub_vec result;
	string line;

	static constexpr size_t ID_OFFSET   = 0;
	static constexpr size_t PROP_OFFSET = 1;

	while ( getline( prm_input, line ) ) {
		const auto   id_itrs   = find_field_itrs( line, ID_OFFSET                                  );
		const auto   prop_itrs = find_field_itrs( line, PROP_OFFSET, 1 + ID_OFFSET, id_itrs.second );
		const auto   id        = make_string_ref        ( id_itrs.first,   id_itrs.second   );
		const auto   props     = parse_double_from_field( prop_itrs.first, prop_itrs.second );
		const size_t id_id     = prm_name_ider.add_name( id );
		if ( result.size() < id_id + 1 ) {
			result.resize( id_id + 1 );
		}
		result[ id_id ] = props;
	}

	return result;
}

/// \brief Parse the cluster item names from the specified string, populating the specified name_ider in the process
doub_vec cath::clust::parse_names(const string      &prm_input,    ///< The string from which to parse the names
                                  id_of_str_bidirnl &prm_name_ider ///< The name_ider to populate from the names data
                                  ) {
	istringstream in_ss{ prm_input };
	return parse_names( in_ss, prm_name_ider );
}

/// \brief Parse the cluster item names from the specified file, populating the specified name_ider in the process
doub_vec cath::clust::parse_names(const path        &prm_input,    ///< The file from which to parse the names
                                  id_of_str_bidirnl &prm_name_ider ///< The name_ider to populate from the names data
                                  ) {
	ifstream in_stream;
	open_ifstream( in_stream, prm_input );
	doub_vec dissims = parse_names( in_stream, prm_name_ider );
	in_stream.close();
	return dissims;
}

/// \brief Parse the cluster item names from the specified istream and return the values and the name_ider
pair<doub_vec, id_of_str_bidirnl> cath::clust::parse_names(istream &prm_input ///< The istream from which to parse the names
                                                           ) {
	id_of_str_bidirnl the_id_of_str_bidirnl;
	auto names = parse_names( prm_input, the_id_of_str_bidirnl );
	return make_pair( std::move( names ), std::move( the_id_of_str_bidirnl ) );
}

/// \brief Parse the cluster item names from the specified string and return the values and the name_ider
pair<doub_vec, id_of_str_bidirnl> cath::clust::parse_names(const string &prm_input ///< The string from which to parse the names
                                                           ) {
	id_of_str_bidirnl the_id_of_str_bidirnl;
	auto names = parse_names( prm_input, the_id_of_str_bidirnl );
	return make_pair( std::move( names ), std::move( the_id_of_str_bidirnl ) );
}

/// \brief Parse the cluster item names from the specified file and return the values and the name_ider
pair<doub_vec, id_of_str_bidirnl> cath::clust::parse_names(const path &prm_input ///< The file from which to parse the names
                                                           ) {
	id_of_str_bidirnl the_id_of_str_bidirnl;
	auto names = parse_names( prm_input, the_id_of_str_bidirnl );
	return make_pair( std::move( names ), std::move( the_id_of_str_bidirnl ) );
}
