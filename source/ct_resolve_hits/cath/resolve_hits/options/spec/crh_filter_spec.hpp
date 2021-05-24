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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_FILTER_SPEC_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_FILTER_SPEC_HPP

#include <boost/algorithm/string/predicate.hpp>
#include <boost/utility/string_ref.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/resolve_hits/hit_score_type.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

#include <regex>

namespace cath::rslv {

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

		/// \brief The (optional) maximum number of queries to process
		size_opt   limit_queries;

		/// \brief The (optional) minimum coverage fraction of an HMM for a hit to be considered
		doub_opt   min_hmm_coverage_frac;

		/// \brief The (optional) minimum coverage fraction of an HMM for a discontinuous hit (/^\dc_{32}$/) to be considered
		doub_opt   min_dc_hmm_coverage_frac;

	public:
		/// \brief The default value for the worst permissible evalue before a hit is ignored
		static constexpr resscr_t DEFAULT_WORST_PERMISSIBLE_EVALUE   = static_cast<resscr_t>( 0.001 );

		/// \brief The default value for the worst permissible bitscore before a hit is ignored
		static constexpr resscr_t DEFAULT_WORST_PERMISSIBLE_BITSCORE = 10.0;

		[[nodiscard]] const resscr_t &  get_worst_permissible_evalue() const;
		[[nodiscard]] const resscr_t &  get_worst_permissible_bitscore() const;
		[[nodiscard]] const resscr_opt &get_worst_permissible_score() const;
		[[nodiscard]] const str_vec &   get_filter_query_ids() const;
		[[nodiscard]] const size_opt &  get_limit_queries() const;
		[[nodiscard]] const doub_opt &  get_min_hmm_coverage_frac() const;
		[[nodiscard]] const doub_opt &  get_min_dc_hmm_coverage_frac() const;

		crh_filter_spec & set_worst_permissible_evalue(const resscr_t &);
		crh_filter_spec & set_worst_permissible_bitscore(const resscr_t &);
		crh_filter_spec & set_worst_permissible_score(const resscr_opt &);
		crh_filter_spec & set_filter_query_ids(const str_vec &);
		crh_filter_spec & set_limit_queries(const size_opt &);
		crh_filter_spec & set_min_hmm_coverage_frac(const doub_opt &);
		crh_filter_spec & set_min_dc_hmm_coverage_frac(const doub_opt &);
	};

	crh_filter_spec make_accept_all_filter_spec();

	/// \brief Return whether the specified score of the specified type passes the specified filter_spec
	inline bool score_passes_filter(const crh_filter_spec &prm_filter_spec, ///< The filter_spec to apply
	                                const double          &prm_score,       ///< The score to test
	                                const hit_score_type  &prm_score_type   ///< The type of score to test
	                                ) {
		switch ( prm_score_type ) {
			case ( hit_score_type::FULL_EVALUE ) : { return  ( prm_score <= prm_filter_spec.get_worst_permissible_evalue  () ); }
			case ( hit_score_type::BITSCORE    ) : { return  ( prm_score >= prm_filter_spec.get_worst_permissible_bitscore() ); }
			case ( hit_score_type::CRH_SCORE   ) : {
				return (
					! prm_filter_spec.get_worst_permissible_score()
					||
					prm_score >= *prm_filter_spec.get_worst_permissible_score() );
			}
		}
		BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_score_type not recognised whilst checking score_passes_filter()"));
	}

	/// \brief Whether the data for the specified query ID should be specified given the specified filter query IDs
	inline bool should_skip_query_id(const str_vec     &prm_filter_query_ids, ///< The list of query IDs to filter (or allow everything if this is empty)
	                                 const std::string &prm_query_id          ///< The query ID string to check
	                                 ) {
		// If there are filter query IDs, then if this amongst them, skip this entry
		return (
			! prm_filter_query_ids.empty()
			&&
			! common::contains( prm_filter_query_ids, prm_query_id )
		);
	}

	/// \brief Whether the data for the specified query ID should be specified given the specified filter query IDs
	inline bool should_skip_query_id(const str_vec           &prm_filter_query_ids, ///< The list of query IDs to filter (or allow everything if this is empty)
	                                 const boost::string_ref &prm_query_id          ///< The query ID string_ref to check
	                                 ) {
		// If there are filter query IDs, then if this amongst them, skip this entry
		return (
			! prm_filter_query_ids.empty()
			&&
			! common::contains( prm_filter_query_ids, prm_query_id )
		);
	}

	/// \brief Return the min-hmm-coverage value to be required of the specified ID under the specified crh_filter_spec
	///        or nullopt if no value is required
	///
	/// \TODO Consider extracting common code from this and cath_score_category_of_id()
	///
	/// \relates crh_filter_spec
	inline doub_opt hmm_coverage_for_match(const crh_filter_spec &prm_filter_spec, ///< The filter_spec to apply
	                                       const std::string     &prm_match_id     ///< The match ID under consideration
	                                       ) {
		static const std::regex  dc_regex        { R"(^dc_\w{32}$)" };
		static const std::string dc_prefix_suffix{ "dc_"            };

		// If a min_dc_hmm_coverage fraction has been specified then that takes precedence if
		// this is a /^dc_\w{32}$/ query ID so check that and return the min_dc_hmm_coverage fraction if so
		const auto &min_dc_hmm_coverage_frac = prm_filter_spec.get_min_dc_hmm_coverage_frac();
		if ( min_dc_hmm_coverage_frac ) {
			if ( prm_match_id.length() == 35 && boost::algorithm::starts_with( prm_match_id, dc_prefix_suffix ) ) {
				if ( regex_search( ::std::cbegin( prm_match_id ), ::std::cend( prm_match_id ), dc_regex ) ) {
					return min_dc_hmm_coverage_frac;
				}
			}
		}

		// Otherwise return any min_hmm_coverage fraction that has been specified (or nullopt if none)
		return prm_filter_spec.get_min_hmm_coverage_frac();
	}

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_OPTIONS_SPEC_CRH_FILTER_SPEC_HPP
