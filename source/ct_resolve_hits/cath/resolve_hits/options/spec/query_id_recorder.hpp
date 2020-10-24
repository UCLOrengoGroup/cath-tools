/// \file
/// \brief The query_id_recorder class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_QUERY_ID_RECORDER_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_QUERY_ID_RECORDER_HPP

#include <boost/utility/string_ref.hpp>

#include "cath/common/algorithm/contains.hpp"

#include <functional>
#include <set>
#include <string>

namespace cath {
	namespace rslv {

		/// \brief A record of the query IDs that have been seen to implement limit_queries
		///
		/// This only needs to be updated when limit_queries has been requested
		///
		/// Use should_skip_query_and_update() rather than calling add_query_id() directly
		///
		/// This uses heterogeneous lookup in an attempt to support lookup with string_ref
		/// without having to construct a string
		class query_id_recorder final {
		private:
			std::set<std::string, std::less<>> seen_query_ids;

		public:
			query_id_recorder() = default;

			size_t size() const;
			bool empty() const;

			bool seen_query_id(const std::string &) const;
			bool seen_query_id(const boost::string_ref &) const;

			query_id_recorder & add_query_id(const std::string &);
			query_id_recorder & add_query_id(const boost::string_ref &);
		};

		/// \brief Return the number of query IDs
		inline size_t query_id_recorder::size() const {
			return seen_query_ids.size();
		}

		/// \brief Return whether this is empty (ie has 0 query IDs)
		inline bool query_id_recorder::empty() const {
			return seen_query_ids.empty();
		}

		/// \brief Return whether the specified query ID has already been seen
		inline bool query_id_recorder::seen_query_id(const std::string &prm_query_id ///< The query ID for which to search
		                                             ) const {
			return common::contains( seen_query_ids, prm_query_id );
		}

		/// \brief Return whether the specified query ID has already been seen
		inline bool query_id_recorder::seen_query_id(const boost::string_ref &prm_query_id ///< The query ID for which to search
		                                             ) const {
			return common::contains( seen_query_ids, prm_query_id );
		}

		/// \brief Add the specified query ID to those that have been seen
		///
		/// Consider using should_skip_query_and_update() rather than calling this method directly
		inline query_id_recorder & query_id_recorder::add_query_id(const std::string &prm_query_id ///< The query ID to add
		                                                           ) {
			seen_query_ids.insert( prm_query_id );
			return *this;
		}

// #define STR_HELPER(x) #x
// #define STR(x) STR_HELPER(x)

// #ifdef __GLIBCPP__
// #pragma message "content of __GLIBCPP__: " STR(__GLIBCPP__)
// #endif

// #ifdef __GLIBCXX__
// #pragma message "content of __GLIBCXX__: " STR(__GLIBCXX__)
// #endif

		/// \brief Add the specified query ID to those that have been seen
		///
		/// Consider using should_skip_query_and_update() rather than calling this method directly
		inline query_id_recorder & query_id_recorder::add_query_id(const boost::string_ref &prm_query_id ///< The query ID to add
		                                                           ) {
			// Find the place where query ID would go
			const auto lower_bound_itr = seen_query_ids.lower_bound(
// workaround absence of hetergeneous lookup in old libstdc++
#ifdef __GLIBCXX__
#if __GLIBCXX__ < 20160801
				static_cast<std::string>(
#endif
#endif
				prm_query_id
#ifdef __GLIBCXX__
#if __GLIBCXX__ < 20160801
				)
#endif
#endif
			);

			// If the query ID should be added, then emplace directly
			if ( lower_bound_itr == common::cend( seen_query_ids ) || *lower_bound_itr != prm_query_id ) {
				seen_query_ids.emplace_hint(
					lower_bound_itr,
					common::cbegin( prm_query_id ),
					common::cend  ( prm_query_id )
				);
			}
			return *this;
		}


	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_QUERY_ID_RECORDER_HPP
