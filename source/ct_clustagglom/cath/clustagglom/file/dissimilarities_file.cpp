/// \file
/// \brief The dissimilarities_file class definitions

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

#include "dissimilarities_file.hpp"

#include <filesystem>
#include <fstream>

#include "cath/clustagglom/link.hpp"
#include "cath/clustagglom/links.hpp"
#include "cath/common/container/id_of_str_bidirnl.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/string/string_parse_tools.hpp"

using namespace ::cath::clust;
using namespace ::cath::common;

using ::std::filesystem::path;
using ::std::ifstream;
using ::std::istream;
using ::std::istringstream;
using ::std::string;

/// \brief Parse the cluster dissimilarities/strengths (ie links) from the specified istream
links cath::clust::parse_dissimilarities(istream           &prm_input,     ///< The istream from which the links should be read
                                         id_of_str_bidirnl &prm_name_ider, ///< The name_ider to populate from the links data
                                         const link_dirn   &prm_link_dirn, ///< Whether the links in the input file represent strengths or dissimilarities
                                         const size_t      &prm_column_idx ///< The index (offset 0) of the column from which the strengths or dissimilarities should be parsed
                                         ) {
	links result;
	string line;

	static constexpr size_t ID1_OFFSET     = 0;
	static constexpr size_t ID2_OFFSET     = 1;

	while ( getline( prm_input, line ) ) {
		const auto     id1_itrs   = find_field_itrs( line, ID1_OFFSET                                      );
		const auto     id2_itrs   = find_field_itrs( line, ID2_OFFSET,     1 + ID1_OFFSET, id1_itrs.second );
		const auto     value_itrs = find_field_itrs( line, prm_column_idx, 1 + ID2_OFFSET, id2_itrs.second );
		const auto     id1        = make_string_ref( id1_itrs.first, id1_itrs.second );
		const auto     id2        = make_string_ref( id2_itrs.first, id2_itrs.second );
		const strength seq_id     = std::is_same_v<strength, float>
		                              ? parse_float_from_field( value_itrs.first, value_itrs.second )
		                              : static_cast<strength>( parse_double_from_field( value_itrs.first, value_itrs.second ) );

		const strength link_val   = ( prm_link_dirn == link_dirn::STRENGTH ) ? -seq_id : seq_id;
		const item_idx id_1_id    = debug_numeric_cast<item_idx>( prm_name_ider.add_name( id1 ) );
		const item_idx id_2_id    = debug_numeric_cast<item_idx>( prm_name_ider.add_name( id2 ) );
		if ( id_1_id != id_2_id ) {
			result.add_link_symmetrically( id_1_id, id_2_id, link_val );
		}
	}

	return result;
}

/// \brief Parse the cluster dissimilarities/strengths (ie links) from the specified string
links cath::clust::parse_dissimilarities(const string      &prm_input,     ///< The string from which the links should be read
                                         id_of_str_bidirnl &prm_name_ider, ///< The name_ider to populate from the links data
                                         const link_dirn   &prm_link_dirn, ///< Whether the links in the input file represent strengths or dissimilarities
                                         const size_t      &prm_column_idx ///< The index (offset 0) of the column from which the strengths or dissimilarities should be parsed
                                         ) {
	istringstream in_ss{ prm_input };
	return parse_dissimilarities( in_ss, prm_name_ider, prm_link_dirn, prm_column_idx );
}

/// \brief Parse the cluster dissimilarities/strengths (ie links) from the specified file
links cath::clust::parse_dissimilarities(const path        &prm_input,     ///< The file from which the links should be read
                                         id_of_str_bidirnl &prm_name_ider, ///< The name_ider to populate from the links data
                                         const link_dirn   &prm_link_dirn, ///< Whether the links in the input file represent strengths or dissimilarities
                                         const size_t      &prm_column_idx ///< The index (offset 0) of the column from which the strengths or dissimilarities should be parsed
                                         ) {
	ifstream in_stream = open_ifstream( prm_input );
	links dissims = parse_dissimilarities( in_stream, prm_name_ider, prm_link_dirn, prm_column_idx );
	in_stream.close();
	return dissims;
}

