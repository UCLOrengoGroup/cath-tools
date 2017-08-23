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

#include "clustagglom/link.hpp"
#include "clustagglom/links.hpp"
#include "common/container/id_of_str_bidirnl.hpp"
#include "common/file/open_fstream.hpp"
#include "common/string/string_parse_tools.hpp"

#include <fstream>

using namespace cath::clust;
using namespace cath::common;

using boost::filesystem::path;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::max;
using std::string;

/// \brief Parse the cluster dissimilarities/strengths (ie links) from the specified istream
links cath::clust::parse_dissimilarities(istream           &arg_input,     ///< The istream from which the links should be read
                                         id_of_str_bidirnl &arg_name_ider, ///< The name_ider to populate from the links data
                                         const link_dirn   &arg_link_dirn  ///< Whether the links in the input file represent strengths or dissimilarities
                                         ) {
	links result;
	string line;

	static constexpr size_t ID1_OFFSET     = 0;
	static constexpr size_t ID2_OFFSET     = 1;
	static constexpr size_t SEQ_ID_OFFSET  = 2;
	static constexpr size_t OVERLAP_OFFSET = 3;

	while ( getline( arg_input, line ) ) {
		const auto     id1_itrs     = find_field_itrs( line, ID1_OFFSET                                                );
		const auto     id2_itrs     = find_field_itrs( line, ID2_OFFSET,     1 + ID1_OFFSET,    id1_itrs.second        );
		const auto     seq_id_itrs  = find_field_itrs( line, SEQ_ID_OFFSET,  1 + ID2_OFFSET,    id2_itrs.second        );
		const auto     overlap_itrs = find_field_itrs( line, OVERLAP_OFFSET, 1 + SEQ_ID_OFFSET, seq_id_itrs.second     );

		if ( parse_float_from_field ( overlap_itrs.first, overlap_itrs.second ) < 80 ) {
			continue;
		}

		const auto     id1          = make_string_ref       ( id1_itrs.first,     id1_itrs.second     );
		const auto     id2          = make_string_ref       ( id2_itrs.first,     id2_itrs.second     );
		const strength seq_id       = std::is_same<strength, float>::value
			? static_cast<strength>( parse_float_from_field ( seq_id_itrs.first, seq_id_itrs.second ) )
			: static_cast<strength>( parse_double_from_field( seq_id_itrs.first, seq_id_itrs.second ) );

		const strength link_val     = ( arg_link_dirn == link_dirn::STRENGTH ) ? -seq_id : seq_id;

		const item_idx id_1_id = debug_numeric_cast<item_idx>( arg_name_ider.add_name( id1 ) );
		const item_idx id_2_id = debug_numeric_cast<item_idx>( arg_name_ider.add_name( id2 ) );
		if ( id_1_id != id_2_id ) {
			result.add_link_symmetrically( id_1_id, id_2_id, link_val );
		}
	}

	return result;
}

/// \brief Parse the cluster dissimilarities/strengths (ie links) from the specified string
links cath::clust::parse_dissimilarities(const string      &arg_input,     ///< The string from which the links should be read
                                         id_of_str_bidirnl &arg_name_ider, ///< The name_ider to populate from the links data
                                         const link_dirn   &arg_link_dirn  ///< Whether the links in the input file represent strengths or dissimilarities
                                         ) {
	istringstream in_ss{ arg_input };
	return parse_dissimilarities( in_ss, arg_name_ider, arg_link_dirn );
}

/// \brief Parse the cluster dissimilarities/strengths (ie links) from the specified file
links cath::clust::parse_dissimilarities(const path        &arg_input,     ///< The file from which the links should be read
                                         id_of_str_bidirnl &arg_name_ider, ///< The name_ider to populate from the links data
                                         const link_dirn   &arg_link_dirn  ///< Whether the links in the input file represent strengths or dissimilarities
                                         ) {
	ifstream in_stream;
	open_ifstream( in_stream, arg_input );
	links dissims = parse_dissimilarities( in_stream, arg_name_ider, arg_link_dirn );
	in_stream.close();
	return dissims;
}

