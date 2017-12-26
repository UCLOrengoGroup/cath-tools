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

#include <boost/log/trivial.hpp>
#include <boost/utility/string_ref.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "cluster/old_cluster_data.hpp"
#include "common/boost_addenda/log/log_to_ostream_guard.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "common/optional/make_optional_if.hpp"
#include "common/string/string_parse_tools.hpp"
#include "seq/seq_seg_run_parser.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace cath;
using namespace cath::clust;
using namespace cath::common;
using namespace cath::seq;

using boost::filesystem::path;
using boost::string_ref;
using std::ifstream;
using std::istream;
using std::istringstream;
using std::ostream;
using std::string;

static constexpr size_t CLUSTER_ID_OFFSET = 0;
static constexpr size_t DOMAIN_ID_OFFSET  = 1;

/// \brief Print any warnings to the specified (optional) ostream arising from the interaction (if any) of the new entry
inline void warn_if_neccessary(const clust_entry_problem &arg_problem,                 ///< The type of problem encountered when reading the new entry
                               const ostream_ref_opt     &arg_ostream_ref_opt,         ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                               const string_ref          &arg_cluster_name,            ///< A string_ref to the name of the cluster
                               const string_ref          &arg_entry_name,              ///< A string_ref to the name of the entry
                               bool                      &arg_warned_duplicate,        ///< Whether a warning about duplicates has already been given
                               const str_opt             &arg_extra_info = boost::none ///< An optional string containing extra information about the problem (if any)
                               ) {
	if ( arg_problem != clust_entry_problem::NONE && arg_ostream_ref_opt ) {
		const log_to_ostream_guard ostream_log_guard{ arg_ostream_ref_opt.get().get() };

		switch ( arg_problem ) {
			case ( clust_entry_problem::REPEAT ) : {
				if ( ! arg_warned_duplicate ) {
					BOOST_LOG_TRIVIAL( warning ) << "Skipping entry "
						<< arg_entry_name   << " (in cluster "
						<< arg_cluster_name << ") because it duplicates a previous entry in the same cluster-membership input data. Will not warn about any further duplicate entries."
						<< ( arg_extra_info ? " - " + *arg_extra_info : string{} );
					arg_warned_duplicate = true;
				}
				break;
			}
			case ( clust_entry_problem::CLASH ) : {
				BOOST_LOG_TRIVIAL( warning ) << "Skipping entry "
					<< arg_entry_name   << " (in cluster "
					<< arg_cluster_name << ") because it clashes with a previous entry in the same cluster-membership input data"
					<< ( arg_extra_info ? " - " + *arg_extra_info : string{} );
				break;
			}
			case ( clust_entry_problem::PARSE_ERROR ) : {
				BOOST_LOG_TRIVIAL( warning ) << "Problem parsing segments from entry "
					<< arg_entry_name   << " (in cluster "
					<< arg_cluster_name << ")"
					<< ( arg_extra_info ? " - " + *arg_extra_info : string{} );
				break;
			}
			case ( clust_entry_problem::NONE ) : {}
		}
	}
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
old_cluster_data cath::clust::parse_old_membership(istream               &arg_istream,           ///< The istream to parse from
                                                   id_of_str_bidirnl     &arg_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &arg_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	using std::to_string;
	bool warned_duplicate = false;
	old_cluster_data result{ arg_id_of_str_bidirnl };
	seq_seg_run_parser segs_parser;
	string line;
	size_t line_ctr = 0;
	while ( getline( arg_istream, line ) ) {
		++line_ctr;
		const auto cluster_id_itrs = find_field_itrs( line, CLUSTER_ID_OFFSET                                               );
		const auto domain_id_itrs  = find_field_itrs( line, DOMAIN_ID_OFFSET, 1 + CLUSTER_ID_OFFSET, cluster_id_itrs.second );

		// Check that the line doesn't have too few or too many fields
		if ( domain_id_itrs.first == domain_id_itrs.second ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				  "Cannot parse cluster membership from fewer than two fields at line number "
				+ to_string( line_ctr )
				+ ". Line is: \""
				+ line
				+ "\""
			));
		}
		if ( find_itr_before_first_non_space( domain_id_itrs.second, common::cend( line ) ) != common::cend( line ) ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				  "Cannot parse cluster membership from more than two fields at line number "
				+ to_string( line_ctr )
				+ ". Line is: \""
				+ line
				+ "\""
			));
		}

		/// \TODO Check there are no more fields (but allow whitespace)
		const auto         slash_index          = make_string_ref( domain_id_itrs ).find_last_of( '/' );
		const bool         has_segs             = ( slash_index != string_ref::npos );
		const auto         pre_split_point_itr  = has_segs ? next( domain_id_itrs.first, debug_numeric_cast<ptrdiff_t>( slash_index     ) )
		                                                   : domain_id_itrs.second;
		const auto cluster_name  = make_string_ref( cluster_id_itrs.first, cluster_id_itrs.second );
		const auto sequence_name = make_string_ref( domain_id_itrs.first,  pre_split_point_itr    );
		const auto entry_name    = make_string_ref( domain_id_itrs.first,  domain_id_itrs.second  );
		try {
			const auto &interaction  = result.add_entry(
				cluster_name,
				sequence_name,
				entry_name,
				make_optional_if_fn(
					has_segs,
					[&] {
						return segs_parser.parse( next( pre_split_point_itr ), domain_id_itrs.second );
					}
				)
			);
			warn_if_neccessary( interaction, arg_ostream_ref_opt, cluster_name, entry_name, warned_duplicate );
		}
		catch (const std::exception &x) {
			warn_if_neccessary(
				clust_entry_problem::PARSE_ERROR,
				arg_ostream_ref_opt,
				cluster_name,
				entry_name,
				warned_duplicate,
				string{ x.what() }
			);
		}
	}
	return result;
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
old_cluster_data cath::clust::parse_old_membership(const string          &arg_input,             ///< The string to parse from
                                                   id_of_str_bidirnl     &arg_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &arg_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	istringstream in_ss{ arg_input };
	return parse_old_membership( in_ss, arg_id_of_str_bidirnl, arg_ostream_ref_opt );
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
old_cluster_data cath::clust::parse_old_membership(const path            &arg_input,             ///< The file to parse from
                                                   id_of_str_bidirnl     &arg_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &arg_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	ifstream in_stream;
	open_ifstream( in_stream, arg_input );
	return parse_old_membership( in_stream, arg_id_of_str_bidirnl, arg_ostream_ref_opt );
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(istream               &arg_istream,           ///< The istream to parse from
                                                   id_of_str_bidirnl     &arg_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &arg_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	using std::to_string;
	bool warned_duplicate = false;
	new_cluster_data result{ arg_id_of_str_bidirnl };
	seq_seg_run_parser segs_parser;
	string line;
	size_t line_ctr = 0;
	while ( getline( arg_istream, line ) ) {
		++line_ctr;
		const auto cluster_id_itrs = find_field_itrs( line, CLUSTER_ID_OFFSET                                               );
		const auto domain_id_itrs  = find_field_itrs( line, DOMAIN_ID_OFFSET, 1 + CLUSTER_ID_OFFSET, cluster_id_itrs.second );

		// Check that the line doesn't have too few or too many fields
		if ( domain_id_itrs.first == domain_id_itrs.second ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				  "Cannot parse cluster membership from fewer than two fields at line number "
				+ to_string( line_ctr )
				+ ". Line is: \""
				+ line
				+ "\""
			));
		}
		if ( find_itr_before_first_non_space( domain_id_itrs.second, common::cend( line ) ) != common::cend( line ) ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				  "Cannot parse cluster membership from more than two fields at line number "
				+ to_string( line_ctr )
				+ ". Line is: \""
				+ line
				+ "\""
			));
		}

		const auto         slash_index          = make_string_ref( domain_id_itrs ).find_last_of( '/' );
		const bool         has_segs             = ( slash_index != string_ref::npos );
		const auto         pre_split_point_itr  = has_segs ? next( domain_id_itrs.first, debug_numeric_cast<ptrdiff_t>( slash_index     ) )
		                                                   : domain_id_itrs.second;
		const auto cluster_name  = make_string_ref( cluster_id_itrs.first, cluster_id_itrs.second );
		const auto sequence_name = make_string_ref( domain_id_itrs.first,  pre_split_point_itr    );
		const auto entry_name    = make_string_ref( domain_id_itrs.first,  domain_id_itrs.second  );
		try {
			const auto problem   = result.add_entry(
				cluster_name,
				sequence_name,
				entry_name,
				make_optional_if_fn(
					has_segs,
					[&] {
						return segs_parser.parse( next( pre_split_point_itr ), domain_id_itrs.second );
					}
				)
			);
			warn_if_neccessary( problem, arg_ostream_ref_opt, cluster_name, entry_name, warned_duplicate );
		}
		catch (const std::exception &x) {
			warn_if_neccessary(
				clust_entry_problem::PARSE_ERROR,
				arg_ostream_ref_opt,
				cluster_name,
				entry_name,
				warned_duplicate,
				string{ x.what() }
			);
		}
	}

	return result;
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(const string          &arg_input,             ///< The string to parse from
                                                   id_of_str_bidirnl     &arg_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &arg_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	istringstream in_ss{ arg_input };
	return parse_new_membership( in_ss, arg_id_of_str_bidirnl, arg_ostream_ref_opt );
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(const path            &arg_input,             ///< The file to parse from
                                                   id_of_str_bidirnl     &arg_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &arg_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	ifstream in_stream;
	open_ifstream( in_stream, arg_input );
	return parse_new_membership( in_stream, arg_id_of_str_bidirnl, arg_ostream_ref_opt );
}
