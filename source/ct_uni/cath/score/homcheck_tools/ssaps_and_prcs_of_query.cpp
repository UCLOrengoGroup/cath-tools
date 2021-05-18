/// \file
/// \brief The ssaps_and_prcs_of_query class definitions

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

#include "ssaps_and_prcs_of_query.hpp"

#include <string>
#include <unordered_map>

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <spdlog/spdlog.h>

#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/range/adaptor/adjacented.hpp"
#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/hash/pair_hash.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/prc_scores_file/prc_scores_entry.hpp"
#include "cath/file/ssap_scores_file/ssap_scores_entry.hpp"
#include "cath/score/homcheck_tools/first_result_if.hpp"
#include "cath/score/homcheck_tools/ssap_and_prc.hpp"
#include "cath/score/homcheck_tools/superfamily_of_domain.hpp"
#include "cath/score/score_classification/rbf_model.hpp"

using namespace ::cath::common;
using namespace ::cath::homcheck;
using namespace ::cath::file;
using namespace ::cath::score;
using namespace ::std;

using ::boost::algorithm::all_of;
using ::boost::range::for_each;

/// \brief Check that the class invariants hold
///
/// \pre All ssap_and_prc_entries have the same query_id else this
///      throws a invalid_argument_exception
void ssaps_and_prcs_of_query::sanity_check() const {
	const bool query_ids_identical = all_of(
		ssap_and_prc_entries | adjacented,
		[] (const pair<const ssap_and_prc &, const ssap_and_prc &> &x) { return x.first.get_query_id() == x.second.get_query_id(); }
	);

	if ( ! query_ids_identical ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Conflicting query IDs in data for ssaps_and_prcs_of_query (first query ID : \""
			+ ssap_and_prc_entries.front().get_query_id()
			+ "\")"
		));
	}
}

/// \brief Ctor from a vector of ssap_and_prc objects
///
/// \pre All ssap_and_prc_entries have the same query_id else this
///      throws a invalid_argument_exception
ssaps_and_prcs_of_query::ssaps_and_prcs_of_query(ssap_and_prc_vec prm_ssap_and_prc_entries ///< The vector of ssap_and_prc objects from which this ssaps_and_prcs_of_query should be constructed
                                                 ) : ssap_and_prc_entries { std::move( prm_ssap_and_prc_entries ) } {
	sanity_check();
}

/// \brief Calculate the SVM scores for all the results using the specified SVM RBF model
void ssaps_and_prcs_of_query::calculate_all_svm_scores(const rbf_model &prm_svm /// \brief The SVM RBF model with which to calculate the scores
                                                       ) {
	for_each(
		ssap_and_prc_entries,
		[&] (ssap_and_prc &x) { x.calculate_svm_score( prm_svm ); }
	);
}

/// \brief Get the whether this ssaps_and_prcs_of_query is empty (ie has no ssap_and_prc objects)
bool ssaps_and_prcs_of_query::empty() const {
	return ssap_and_prc_entries.empty();
}

/// \brief Get the size of this ssaps_and_prcs_of_query (ie the number of ssap_and_prc objects)
size_t ssaps_and_prcs_of_query::size() const {
	return ssap_and_prc_entries.size();
}

/// \brief Get the ssap_and_prc object at the specified index
const ssap_and_prc & ssaps_and_prcs_of_query::operator[](const size_t &prm_index ///< The index of the ssap_and_prc object to access
                                                         ) const {
	return ssap_and_prc_entries[ prm_index ];
}

/// \brief Standard begin() function to make this into a const range
auto ssaps_and_prcs_of_query::begin() const -> const_iterator {
	return cbegin( ssap_and_prc_entries );
}

/// \brief Standard end() function to make this into a const range
auto ssaps_and_prcs_of_query::end() const -> const_iterator{
	return cend( ssap_and_prc_entries );
}

/// \brief Calculate the SVM scores for all the results in a copy of the specified ssaps_and_prcs_of_query and return that copy
ssaps_and_prcs_of_query cath::homcheck::calculate_all_svm_scores_copy(ssaps_and_prcs_of_query  prm_ssaps_and_prcs, ///< The ssaps_and_prcs_of_query from which a copy should be taken, updated and returned
                                                                      const rbf_model         &prm_svm             ///< The SVM RBF model with which to calculate the scores
                                                                      ) {
	prm_ssaps_and_prcs.calculate_all_svm_scores( prm_svm );
	return prm_ssaps_and_prcs;
}

/// \brief Get the query ID from the specified SSAP and PRC results
///
/// \relates ssaps_and_prcs_of_query
const string & cath::homcheck::get_query_id(const ssaps_and_prcs_of_query &prm_ssaps_and_prcs /// The SSAP and PRC results to query
                                            ) {
	if ( prm_ssaps_and_prcs.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Can't get query of empty ssaps_and_prcs_of_query object"));
	}
	return front( prm_ssaps_and_prcs ).get_query_id();
}

/// \brief Return the best by SVM hit to a domain in CATH of the specified results
///
/// \relates ssaps_and_prcs_of_query
ssap_and_prc_cref_opt cath::homcheck::best_svm_assignable(const ssaps_and_prcs_of_query &prm_ssaps_and_prcs,        ///< The SSAP and PRC results for the query domain
                                                          const superfamily_of_domain   &prm_superfamily_of_domain, ///< The superfamily_of_domain for finding which matches are assigned
                                                          const double                  &prm_minimum_ssap_overlap   ///< TODOCUMENT
                                                          ) {
	return first_result_if(
		prm_ssaps_and_prcs,
		[&] (const ssap_and_prc &x, const ssap_and_prc &y) {
			// Reverse inequality to put the highest magic_function values to the start
			return get_svm_score( x ) > get_svm_score( y );
		},
		[&] (const ssap_and_prc &x) {
			return (
				get_svm_score( x ) >= 3.47554072714
				&&
				get_ssap_overlap_pc( x ) >= prm_minimum_ssap_overlap
				&&
				prm_superfamily_of_domain.has_superfamily_of_domain( x.get_match_id() )
			);
		}
	);
}

/// \brief Return the best by magic-function hit to a domain in CATH of the specified results
///
/// \relates ssaps_and_prcs_of_query
ssap_and_prc_cref_opt cath::homcheck::best_magic_function_assignable(const ssaps_and_prcs_of_query &prm_ssaps_and_prcs,        ///< The SSAP and PRC results for the query domain
                                                                     const superfamily_of_domain   &prm_superfamily_of_domain, ///< The superfamily_of_domain for finding which matches are assigned
                                                                     const double                  &prm_minimum_ssap_overlap   ///< TODOCUMENT
                                                                     ) {
	return first_result_if(
		prm_ssaps_and_prcs,
		[] (const ssap_and_prc &x, const ssap_and_prc &y) {
			// Reverse inequality to put the highest magic_function values to the start
			return x.get_magic_function_score() > y.get_magic_function_score();
		},
		[&] (const ssap_and_prc &x) {
			return (
				x.get_magic_function_score() >= 79.7779328254
				&&
				get_ssap_overlap_pc( x ) >= prm_minimum_ssap_overlap
				&&
				prm_superfamily_of_domain.has_superfamily_of_domain( x.get_match_id() )
			);
		}
	);
}

/// \brief The best fold-level match to a domain in CATH from the specified SSAP results
///
/// This uses Christine's criteria for a fold level match (SSAP score >= 70.0 && SSAP overlap >= 60.0)
///
/// \relates ssaps_and_prcs_of_query
ssap_scores_entry_cref_opt cath::homcheck::best_fold_level_match(const ssap_scores_entry_vec &prm_ssaps,                ///< The SSAP results for the query domain
                                                                 const superfamily_of_domain &prm_superfamily_of_domain ///< The superfamily_of_domain for finding which matches are assigned
                                                                 ) {
	return first_result_if(
		prm_ssaps,
		[] (const ssap_scores_entry &x, const ssap_scores_entry &y) {
			// Reverse inequality to put the highest SSAP scores to the start
			return x.get_ssap_score() > y.get_ssap_score();
		},
		[&] (const ssap_scores_entry &x) {
			return (
				x.get_ssap_score() >= 70.0
				&&
				x.get_overlap_pc() >= 60.0
				&&
				prm_superfamily_of_domain.has_superfamily_of_domain( x.get_name_2() )
			);
		}
	);
}

/// \brief Build a ssaps_and_prcs_of_query from the specified ssap_scores_entries and prc_scores_entries
///
/// \relates ssaps_and_prcs_of_query
ssaps_and_prcs_of_query cath::homcheck::make_ssaps_and_prcs_of_query(const ssap_scores_entry_vec &prm_ssaps, ///< The SSAPs from which the ssaps_and_prcs_of_query should be built
                                                                     const prc_scores_entry_vec  &prm_prcs   ///< The PRCs from which the ssaps_and_prcs_of_query should be built
                                                                     ) {
	// Sanity check the inputs - step 1: check SSAP results have identical name_1s
	const bool ssap_query_ids_identical = all_of(
		prm_ssaps | adjacented,
		[] (const pair<const ssap_scores_entry &, const ssap_scores_entry &> &x) { return x.first.get_name_1() == x.second.get_name_1(); }
	);
	if ( ! ssap_query_ids_identical ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot construct ssaps_and_prcs_of_query from data because the SSAP results have inconsistent query IDs (first is \""
			+ prm_ssaps.front().get_name_1()
			+ "\")"
		));
	}

	// Sanity check the inputs - step 2: check PRC results have identical name_1s
	const bool prc_query_ids_identical = all_of(
		prm_prcs | adjacented,
		[] (const pair<const prc_scores_entry &, const prc_scores_entry &> &x) { return x.first.get_name_1() == x.second.get_name_1(); }
	);
	if ( ! prc_query_ids_identical ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct ssaps_and_prcs_of_query from data because the PRC results have inconsistent query IDs (first is \""
			+ prm_prcs.front().get_name_1()
			+ "\")"
		));
	}

	// Sanity check the inputs - step 3: check the SSAP name_1 matches the PRC name_1
	if ( ! prm_ssaps.empty() && ! prm_prcs.empty() ) {
		const string &ssaps_name_1 = prm_ssaps.front().get_name_1();
		const string &prcs_name_1  =  prm_prcs.front().get_name_1();
		if ( ssaps_name_1 != prcs_name_1 ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Cannot construct ssaps_and_prcs_of_query from data because the SSAP results query ID \""
				+ ssaps_name_1
				+ "\" does not match the PRC results query ID \""
				+ prcs_name_1
				+ "\""
			));
		}
	}

	const auto prc_index_of_ids = transform_build<unordered_map<str_str_pair, size_t, pair_hash> >(
		indices( prm_prcs.size() ),
		[&] (const size_t &x) {
			const auto &the_prc = prm_prcs[ x ];
			return make_pair( make_pair( the_prc.get_name_1(), the_prc.get_name_2() ), x );
		}
	);

	vector<ssap_and_prc> ssap_and_prc_entries;
	for (const ssap_scores_entry &the_ssap : prm_ssaps) {
		const string &query_id = the_ssap.get_name_1();
		const string &match_id = the_ssap.get_name_2();

		const auto prc_index_itr = prc_index_of_ids.find( make_pair( query_id, match_id ) );
		if ( prc_index_itr != cend( prc_index_of_ids ) ) {
			ssap_and_prc_entries.emplace_back(
				the_ssap,
				prm_prcs[ prc_index_itr->second ]
			);
		}
	}

	const auto num_ssaps = prm_ssaps.size();
	const auto num_prcs  = prm_prcs.size();
	const auto num_comb  = ssap_and_prc_entries.size();
	const string query_id_str        = ( ! ssap_and_prc_entries.empty() ) ? ( " (query: " + get_query_id( ssaps_and_prcs_of_query{ ssap_and_prc_entries } ) + ")" ) : "";
	const string unmatched_ssaps_str = ( num_ssaps > num_comb           ) ? ( std::to_string( num_ssaps - num_comb ) + " unmatched SSAP results from " + std::to_string( num_ssaps ) ) : "";
	const string unmatched_prcs_str  = ( num_prcs  > num_comb           ) ? ( std::to_string( num_prcs  - num_comb ) + " unmatched PRC results from "  + std::to_string( num_prcs  ) ) : "";
	const string conjuction_str      = ( unmatched_ssaps_str.empty() || unmatched_prcs_str.empty() ) ? "" : " and ";
	if ( ! unmatched_ssaps_str.empty() || ! unmatched_prcs_str.empty() ) {
		::spdlog::warn( "After parsing {} ssaps_and_prcs_of_query{}, was left with {}{}{}",
		                num_comb,
		                query_id_str,
		                unmatched_ssaps_str,
		                conjuction_str,
		                unmatched_prcs_str );
	}
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return ssaps_and_prcs_of_query{ ssap_and_prc_entries };
}
