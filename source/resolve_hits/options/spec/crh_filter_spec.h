/// \file
/// \brief The crh_filter_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_FILTER_SPEC_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_FILTER_SPEC_H

#include <boost/optional.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/algorithm/contains.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "resolve_hits/hit_score_type.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

namespace cath {
	namespace rslv {

		/// \brief Specify what filtering should be done of incoming hits
		///
		/// Note that filter_query_ids filters out whole result sets based
		/// on query_id whereas the worst_permissible attributes filter out individual hits
		class crh_filter_spec final {
		private:
			/// \brief The worst permissible evalue before a hit is ignored
			resscr_t   worst_permissible_evalue   = DEFAULT_WORST_PERMISSIBLE_EVALUE;

			/// \brief The worst permissible bitscore before a hit is ignored
			resscr_t   worst_permissible_bitscore = DEFAULT_WORST_PERMISSIBLE_BITSCORE;

			/// \brief The worst permissible cath-resolve-hits score before a hit is ignored
			resscr_opt worst_permissible_score;

			/// \brief The query IDs on which to filter the input, if any are present
			str_vec    filter_query_ids;

		public:
			/// \brief The default value for the worst permissible evalue before a hit is ignored
			static constexpr resscr_t DEFAULT_WORST_PERMISSIBLE_EVALUE   = static_cast<resscr_t>( 0.001 );

			/// \brief The default value for the worst permissible bitscore before a hit is ignored
			static constexpr resscr_t DEFAULT_WORST_PERMISSIBLE_BITSCORE = 10.0;

			const resscr_t & get_worst_permissible_evalue() const;
			const resscr_t & get_worst_permissible_bitscore() const;
			const resscr_opt & get_worst_permissible_score() const;
			const str_vec & get_filter_query_ids() const;

			crh_filter_spec & set_worst_permissible_evalue(const resscr_t &);
			crh_filter_spec & set_worst_permissible_bitscore(const resscr_t &);
			crh_filter_spec & set_worst_permissible_score(const resscr_opt &);
			crh_filter_spec & set_filter_query_ids(const str_vec &);
		};

		crh_filter_spec make_accept_all_filter_spec();

		/// \brief Return whether the specified score of the specified type passes the specified filter_spec
		inline bool score_passes_filter(const crh_filter_spec &arg_filter_spec, ///< The filter_spec to apply
		                                const double          &arg_score,       ///< The score to test
		                                const hit_score_type  &arg_score_type   ///< The type of score to test
		                                ) {
			switch ( arg_score_type ) {
				case ( hit_score_type::FULL_EVALUE ) : { return  ( arg_score <= arg_filter_spec.get_worst_permissible_evalue  () ); }
				case ( hit_score_type::BITSCORE    ) : { return  ( arg_score >= arg_filter_spec.get_worst_permissible_bitscore() ); }
				case ( hit_score_type::CRH_SCORE   ) : {
					return (
						! arg_filter_spec.get_worst_permissible_score()
						||
						arg_score >= *arg_filter_spec.get_worst_permissible_score() );
				}
				default : {
					BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_score_type not recognised whilst checking score_passes_filter()"));
					return true; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
				}
			}
		}

		/// \brief Whether the data for the specified query ID should be specified given the specified filter query IDs
		inline bool should_skip_query_id(const str_vec     &arg_filter_query_ids, ///< The list of query IDs to filter (or allow everything if this is empty)
		                                 const std::string &arg_query_id          ///< The query ID string to check
		                                 ) {
			// If there are filter query IDs, then if this amongst them, skip this entry
			return (
				! arg_filter_query_ids.empty()
				&&
				! common::contains( arg_filter_query_ids, arg_query_id )
			);
		}

		/// \brief Whether the data for the specified query ID should be specified given the specified filter query IDs
		inline bool should_skip_query_id(const str_vec           &arg_filter_query_ids, ///< The list of query IDs to filter (or allow everything if this is empty)
		                                 const boost::string_ref &arg_query_id          ///< The query ID string_ref to check
		                                 ) {
			// If there are filter query IDs, then if this amongst them, skip this entry
			return (
				! arg_filter_query_ids.empty()
				&&
				! common::contains( arg_filter_query_ids, arg_query_id )
			);
		}

		/// \brief Whether the data for the specified query ID should be specified given the specified filter query IDs
		///
		/// \relates crh_filter_spec
		inline bool should_skip_query_id(const crh_filter_spec &arg_filter_spec, ///< The filter_spec to apply
		                                 const std::string     &arg_query_id     ///< The query ID string to check
		                                 ) {
			return should_skip_query_id( arg_filter_spec.get_filter_query_ids(), arg_query_id );
		}

		/// \brief Whether the data for the specified query ID should be specified given the specified filter query IDs
		///
		/// \relates crh_filter_spec
		inline bool should_skip_query_id(const crh_filter_spec   &arg_filter_spec, ///< The filter_spec to apply
		                                 const boost::string_ref &arg_query_id     ///< The query ID string_ref to check
		                                 ) {
			return should_skip_query_id( arg_filter_spec.get_filter_query_ids(), arg_query_id );
		}
		

	} // namespace rslv
} // namespace cath

#endif
