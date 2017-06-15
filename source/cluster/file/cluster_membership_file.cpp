/// \file
/// \brief The cluster_membership_file class definitions

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

#include "cluster_membership_file.hpp"

#include <boost/utility/string_ref.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/file/open_fstream.hpp"
#include "common/optional/make_optional_if.hpp"
#include "common/string/string_parse_tools.hpp"
#include "seq/seq_seg_run_parser.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace cath::clust;
using namespace cath::common;
using namespace cath::seq;

using boost::filesystem::path;
using boost::string_ref;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::string;

/// \brief Parse the data for old "from" clusters from a cluster membership istream
int cath::clust::parse_old_membership(istream &arg_istream ///< The istream to parse from
                                      ) {
	seq_seg_run_parser segs_parser;
	string line;
	while ( getline( arg_istream, line ) ) {
		static constexpr size_t CLUSTER_ID_OFFSET = 0;
		static constexpr size_t DOMAIN_ID_OFFSET  = 1;
		const auto cluster_id_itrs = find_field_itrs( line, CLUSTER_ID_OFFSET                                               );
		const auto domain_id_itrs  = find_field_itrs( line, DOMAIN_ID_OFFSET, 1 + CLUSTER_ID_OFFSET, cluster_id_itrs.second );
		// TODO Check there are no more fields (but allow whitespace)

		const cluster_id_t cluster_id           = parse_uint_from_field( cluster_id_itrs.first, cluster_id_itrs.second );
		const auto         domain_id_str_ref    = make_string_ref( domain_id_itrs );
		const auto         slash_index          = domain_id_str_ref.find_last_of( '/' );
		const bool         has_segs             = ( slash_index != string_ref::npos );
		const auto         pre_split_point_itr  = has_segs ? next( domain_id_itrs.first, debug_numeric_cast<ptrdiff_t>( slash_index     ) )
		                                                   : domain_id_itrs.second;
		const auto         seq_id_str_ref       = make_string_ref( domain_id_itrs.first, pre_split_point_itr );
		const auto         opt_segs             = make_optional_if_fn(
			has_segs,
			[&] () {
				return segs_parser.parse( next( pre_split_point_itr ), domain_id_itrs.second );
			}
		);

		ignore_unused( seq_id_str_ref, cluster_id );
	}
	return 0;
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
int cath::clust::parse_old_membership(const string &arg_input ///< The string to parse from
                                      ) {
	istringstream in_ss{ arg_input };
	return parse_old_membership( in_ss );
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
int cath::clust::parse_old_membership(const path &arg_input ///< The file to parse from
                                      ) {
	ifstream in_stream;
	open_ifstream( in_stream, arg_input );
	return parse_old_membership( in_stream );
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(istream      &arg_istream,     ///< The istream to parse from
                                                   id_of_string &arg_id_of_string ///< The id_of_string to use to map from sequences names to IDs 
                                                   ) {
	new_cluster_data result{ arg_id_of_string };
	seq_seg_run_parser segs_parser;
	string line;
	while ( getline( arg_istream, line ) ) {
		static constexpr size_t CLUSTER_ID_OFFSET = 0;
		static constexpr size_t DOMAIN_ID_OFFSET  = 1;
		const auto cluster_id_itrs = find_field_itrs( line, CLUSTER_ID_OFFSET                                               );
		const auto domain_id_itrs  = find_field_itrs( line, DOMAIN_ID_OFFSET, 1 + CLUSTER_ID_OFFSET, cluster_id_itrs.second );
		// TODO Check there are no more fields (but allow whitespace)

		const auto         slash_index          = make_string_ref( domain_id_itrs ).find_last_of( '/' );
		const bool         has_segs             = ( slash_index != string_ref::npos );
		const auto         pre_split_point_itr  = has_segs ? next( domain_id_itrs.first, debug_numeric_cast<ptrdiff_t>( slash_index     ) )
		                                                   : domain_id_itrs.second;
		result.add_entry(
			make_string_ref( cluster_id_itrs.first, cluster_id_itrs.second ),
			make_string_ref( domain_id_itrs.first,  pre_split_point_itr    ),
			make_optional_if_fn(
				has_segs,
				[&] () {
					return segs_parser.parse( next( pre_split_point_itr ), domain_id_itrs.second );
				}
			)
		);
	}

	return result;
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(const string &arg_input,       ///< The string to parse from
                                                   id_of_string &arg_id_of_string ///< The id_of_string to use to map from sequences names to IDs
                                                   ) {
	istringstream in_ss{ arg_input };
	return parse_new_membership( in_ss, arg_id_of_string );

}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(const path   &arg_input,       ///< The file to parse from
                                                   id_of_string &arg_id_of_string ///< The id_of_string to use to map from sequences names to IDs
                                                   ) {
	ifstream in_stream;
	open_ifstream( in_stream, arg_input );
	return parse_new_membership( in_stream, arg_id_of_string );
}

