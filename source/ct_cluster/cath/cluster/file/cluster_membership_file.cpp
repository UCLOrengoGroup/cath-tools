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

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>

#include <boost/utility/string_ref.hpp>

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include "cath/cluster/cluster_type_aliases.hpp"
#include "cath/cluster/old_cluster_data.hpp"
#include "cath/common/boost_addenda/log/log_to_ostream_guard.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/optional/make_optional_if.hpp"
#include "cath/common/string/string_parse_tools.hpp"
#include "cath/seq/seq_seg_run_parser.hpp"

using namespace ::cath;
using namespace ::cath::clust;
using namespace ::cath::common;
using namespace ::cath::seq;

using ::boost::string_ref;
using ::std::filesystem::path;
using ::std::ifstream;
using ::std::istream;
using ::std::istringstream;
using ::std::make_optional;
using ::std::nullopt;
using ::std::ostream;
using ::std::string;

static constexpr size_t CLUSTER_ID_OFFSET = 0;
static constexpr size_t DOMAIN_ID_OFFSET  = 1;

/// \brief Print any warnings to the specified (optional) ostream arising from the interaction (if any) of the new entry
static inline void warn_if_necessary(const clust_entry_problem &prm_problem,                    ///< The type of problem encountered when reading the new entry
                                     const ostream_ref_opt     &prm_ostream_ref_opt,            ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                     const string_ref          &prm_cluster_name,               ///< A string_ref to the name of the cluster
                                     const string_ref          &prm_entry_name,                 ///< A string_ref to the name of the entry
                                     bool                      &prm_warned_duplicate,           ///< Whether a warning about duplicates has already been given
                                     const str_opt             &prm_extra_info = ::std::nullopt ///< An optional string containing extra information about the problem (if any)
                                     ) {
	if ( prm_problem != clust_entry_problem::NONE && prm_ostream_ref_opt ) {
		const log_to_ostream_guard ostream_log_guard{ prm_ostream_ref_opt->get() };

		switch ( prm_problem ) {
			case ( clust_entry_problem::REPEAT ) : {
				if ( ! prm_warned_duplicate ) {
					::spdlog::warn(
					  "Skipping entry {} (in cluster {}) because it duplicates a previous entry in the same "
					  "cluster-membership input data. Will not warn about any further duplicate entries.{}",
					  prm_entry_name,
					  prm_cluster_name,
					  ( prm_extra_info ? " - " + *prm_extra_info : string{} ) );
					prm_warned_duplicate = true;
				}
				break;
			}
			case ( clust_entry_problem::CLASH ) : {
				::spdlog::warn( "Skipping entry {} (in cluster {}) because it clashes with a previous entry in the "
				                "same cluster-membership input data{}",
				                prm_entry_name,
				                prm_cluster_name,
				                ( prm_extra_info ? " - " + *prm_extra_info : string{} ) );
				break;
			}
			case ( clust_entry_problem::PARSE_ERROR ) : {
				::spdlog::warn( "Problem parsing segments from entry {} (in cluster {}){}",
				                prm_entry_name,
				                prm_cluster_name,
				                ( prm_extra_info ? " - " + *prm_extra_info : string{} ) );
				break;
			}
			case ( clust_entry_problem::NONE ) : {}
		}
	}
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
old_cluster_data cath::clust::parse_old_membership(istream               &prm_istream,           ///< The istream to parse from
                                                   id_of_str_bidirnl     &prm_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &prm_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	using ::std::to_string;
	bool warned_duplicate = false;
	old_cluster_data result{ prm_id_of_str_bidirnl };
	seq_seg_run_parser segs_parser;
	string line;
	size_t line_ctr = 0;
	while ( getline( prm_istream, line ) ) {
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
		if ( find_itr_before_first_non_space( domain_id_itrs.second, cend( line ) ) != cend( line ) ) {
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
				if_then_optional(
					has_segs,
					segs_parser.parse( next( pre_split_point_itr ), domain_id_itrs.second )
				)
			);
			warn_if_necessary( interaction, prm_ostream_ref_opt, cluster_name, entry_name, warned_duplicate );
		}
		catch (const std::exception &x) {
			warn_if_necessary(
				clust_entry_problem::PARSE_ERROR,
				prm_ostream_ref_opt,
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
old_cluster_data cath::clust::parse_old_membership(const string          &prm_input,             ///< The string to parse from
                                                   id_of_str_bidirnl     &prm_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &prm_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	istringstream in_ss{ prm_input };
	return parse_old_membership( in_ss, prm_id_of_str_bidirnl, prm_ostream_ref_opt );
}

/// \brief Parse the data for old "from" clusters from a cluster membership istream
old_cluster_data cath::clust::parse_old_membership(const path            &prm_input,             ///< The file to parse from
                                                   id_of_str_bidirnl     &prm_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &prm_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	ifstream in_stream = open_ifstream( prm_input );
	return parse_old_membership( in_stream, prm_id_of_str_bidirnl, prm_ostream_ref_opt );
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(istream               &prm_istream,           ///< The istream to parse from
                                                   id_of_str_bidirnl     &prm_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &prm_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	using ::std::to_string;
	bool warned_duplicate = false;
	new_cluster_data result{ prm_id_of_str_bidirnl };
	seq_seg_run_parser segs_parser;
	string line;
	size_t line_ctr = 0;
	while ( getline( prm_istream, line ) ) {
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
		if ( find_itr_before_first_non_space( domain_id_itrs.second, cend( line ) ) != cend( line ) ) {
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
				if_then_optional(
					has_segs,
					segs_parser.parse( next( pre_split_point_itr ), domain_id_itrs.second )
				)
			);
			warn_if_necessary( problem, prm_ostream_ref_opt, cluster_name, entry_name, warned_duplicate );
		}
		catch (const std::exception &x) {
			warn_if_necessary(
				clust_entry_problem::PARSE_ERROR,
				prm_ostream_ref_opt,
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
new_cluster_data cath::clust::parse_new_membership(const string          &prm_input,             ///< The string to parse from
                                                   id_of_str_bidirnl     &prm_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &prm_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	istringstream in_ss{ prm_input };
	return parse_new_membership( in_ss, prm_id_of_str_bidirnl, prm_ostream_ref_opt );
}

/// \brief Parse the data for new "to" clusters from a cluster membership istream
new_cluster_data cath::clust::parse_new_membership(const path            &prm_input,             ///< The file to parse from
                                                   id_of_str_bidirnl     &prm_id_of_str_bidirnl, ///< The id_of_str_bidirnl to use to map from sequences names to IDs
                                                   const ostream_ref_opt &prm_ostream_ref_opt    ///< An optional ostream ref to which warnings about parsing (eg duplicates/clashes) can be written
                                                   ) {
	ifstream in_stream = open_ifstream( prm_input );
	return parse_new_membership( in_stream, prm_id_of_str_bidirnl, prm_ostream_ref_opt );
}
