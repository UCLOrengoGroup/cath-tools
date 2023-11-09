/// \file
/// \brief The should_skip_query class header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_SHOULD_SKIP_QUERY_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_SHOULD_SKIP_QUERY_HPP

#include "cath/resolve_hits/options/spec/crh_filter_spec.hpp"
#include "cath/resolve_hits/options/spec/query_id_recorder.hpp"

namespace cath::rslv {

	/// \brief Return whether the specified query ID should be skipped based limiting the number of queries
	///        to the specified (optional) number given the specified record of previously seen query IDs
	inline bool should_skip_query_by_num(const size_opt          &prm_limit_queries, ///< The (optional) number of queries to which processing should be limited
	                                     const std::string       &prm_query_id,      ///< The query ID string to check
	                                     const query_id_recorder &prm_seen_queries   ///< The query IDs that have already been seen
	                                     ) {
		return (
			prm_limit_queries
			&&
			prm_seen_queries.size() >= *prm_limit_queries
			&&
			! prm_seen_queries.seen_query_id( prm_query_id )
		);
	}

	/// \brief Return whether the specified query ID should be skipped based limiting the number of queries
	///        to the specified (optional) number given the specified record of previously seen query IDs
	inline bool should_skip_query_by_num(const size_opt           &prm_limit_queries, ///< The (optional) number of queries to which processing should be limited
	                                     const ::std::string_view &prm_query_id,      ///< The query ID string_view to check
	                                     const query_id_recorder  &prm_seen_queries   ///< The query IDs that have already been seen
	                                     ) {
		return (
			prm_limit_queries
			&&
			prm_seen_queries.size() >= *prm_limit_queries
			&&
			! prm_seen_queries.seen_query_id( prm_query_id )
		);
	}

	/// \brief Whether the data for the specified query ID should be skipped given the specified filter query IDs
	///
	/// \relates crh_filter_spec
	inline bool should_skip_query(const crh_filter_spec   &prm_filter_spec, ///< The filter_spec to apply
	                              const std::string       &prm_query_id,    ///< The query ID string to check
	                              const query_id_recorder &prm_seen_queries ///< The query IDs that have already been seen
	                              ) {
		return
			should_skip_query_id    ( prm_filter_spec.get_filter_query_ids(), prm_query_id                   )
			||
			should_skip_query_by_num( prm_filter_spec.get_limit_queries(),    prm_query_id, prm_seen_queries );
	}

	/// \brief Whether the data for the specified query ID should be skipped given the specified filter query IDs
	///
	/// \relates crh_filter_spec
	inline bool should_skip_query(const crh_filter_spec    &prm_filter_spec, ///< The filter_spec to apply
	                              const ::std::string_view &prm_query_id,    ///< The query ID string_view to check
	                              const query_id_recorder  &prm_seen_queries ///< The query IDs that have already been seen
	                              ) {
		return
			should_skip_query_id    ( prm_filter_spec.get_filter_query_ids(), prm_query_id                   )
			||
			should_skip_query_by_num( prm_filter_spec.get_limit_queries(),    prm_query_id, prm_seen_queries );
	}

	/// \brief Update the specified record of seen query IDs with the specified query ID if
	///        that's required to implement the specified crh_filter_spec
	///
	/// \relates crh_filter_spec
	inline void update_seen_queries_if_relevant(const crh_filter_spec &prm_filter_spec, ///< The filter_spec to apply
	                                            const std::string     &prm_query_id,    ///< The query ID string to check
	                                            query_id_recorder     &prm_seen_queries ///< The query IDs that have already been seen
	                                            ) {
		if ( prm_filter_spec.get_limit_queries() ) {
			prm_seen_queries.add_query_id( prm_query_id );
		};
	}

	/// \brief Update the specified record of seen query IDs with the specified query ID if
	///        that's required to implement the specified crh_filter_spec
	///
	/// \relates crh_filter_spec
	inline void update_seen_queries_if_relevant(const crh_filter_spec    &prm_filter_spec, ///< The filter_spec to apply
	                                            const ::std::string_view &prm_query_id,    ///< The query ID string_view to check
	                                            query_id_recorder        &prm_seen_queries ///< The query IDs that have already been seen
	                                            ) {
		if ( prm_filter_spec.get_limit_queries() ) {
			prm_seen_queries.add_query_id( prm_query_id );
		};
	}

	/// \brief Return whether the specified query ID should be skipped given the specified crh_filter_spec
	///        and, if appropriate, update the specified list of seen query IDs
	///
	/// \relates crh_filter_spec
	inline bool should_skip_query_and_update(const crh_filter_spec &prm_filter_spec, ///< The filter_spec to apply
	                                         const std::string     &prm_query_id,    ///< The query ID string to check
	                                         query_id_recorder     &prm_seen_queries ///< The query IDs that have already been seen
	                                         ) {
		const bool should_skip = should_skip_query( prm_filter_spec, prm_query_id, prm_seen_queries );
		if ( ! should_skip ) {
			update_seen_queries_if_relevant( prm_filter_spec, prm_query_id, prm_seen_queries );
		}
		return should_skip;
	}

	/// \brief Return whether the specified query ID should be skipped given the specified crh_filter_spec
	///        and, if appropriate, update the specified list of seen query IDs
	///
	/// \relates crh_filter_spec
	inline bool should_skip_query_and_update(const crh_filter_spec    &prm_filter_spec, ///< The filter_spec to apply
	                                         const ::std::string_view &prm_query_id,    ///< The query ID string_view to check
	                                         query_id_recorder        &prm_seen_queries ///< The query IDs that have already been seen
	                                         ) {
		const bool should_skip = should_skip_query( prm_filter_spec, prm_query_id, prm_seen_queries );
		if ( ! should_skip ) {
			update_seen_queries_if_relevant( prm_filter_spec, prm_query_id, prm_seen_queries );
		}
		return should_skip;
	}

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_SHOULD_SKIP_QUERY_HPP
