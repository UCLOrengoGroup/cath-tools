/// \file
/// \brief The ssaps_and_prcs_of_query class definitions

#include "ssaps_and_prcs_of_query.h"

#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/log/trivial.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/irange.hpp>

#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/range/adaptor/adjacented.h"
#include "common/size_t_literal.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "file/prc_scores_file/prc_scores_entry.h"
#include "file/ssap_scores_file/ssap_scores_entry.h"
#include "score/homcheck_tools/ssap_and_prc.h"

//#include <iostream> // ***** TEMPORARY? *****
#include <string>
#include <unordered_map>

using namespace cath::common;
using namespace cath::homcheck;
using namespace cath::file;
using namespace std;

using boost::algorithm::all_of;
using boost::irange;
using boost::make_optional;
using boost::none;
using boost::optional;
using boost::range::max_element;

/// \brief Hashing for pairs of elements
///
/// This is a general class and should be moved out
/// of this file if it can also be helpful elsewhere
struct pair_hash {
public:
	/// \brief Generic function operator for hashing a pair by
	///        combining the results of a standard hash of each value
	template <typename T, typename U>
	size_t operator()(const pair<T, U> &arg_pair ///< The
	                  ) const {
		return hash<T>()( arg_pair.first ) ^ std::hash<U>()( arg_pair.second );
	}
};

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
ssaps_and_prcs_of_query::ssaps_and_prcs_of_query(const ssap_and_prc_vec &arg_ssap_and_prc_entries ///< The vector of ssap_and_prc objects from which this ssaps_and_prcs_of_query should be constructed
                                                 ) : ssap_and_prc_entries ( arg_ssap_and_prc_entries ) {
	sanity_check();
}

/// \brief Get the size of this ssaps_and_prcs_of_query (ie the number of ssap_and_prc objects)
size_t ssaps_and_prcs_of_query::size() const {
	return ssap_and_prc_entries.size();
}

/// \brief Get the ssap_and_prc object at the specified index
const ssap_and_prc & ssaps_and_prcs_of_query::operator[](const size_t &arg_index ///< The index of the ssap_and_prc object to access
                                                         ) const {
	return ssap_and_prc_entries[ arg_index ];
}

/// \brief Standard begin() function to make this into a const range
auto ssaps_and_prcs_of_query::begin() const -> const_iterator {
	return std::cbegin( ssap_and_prc_entries );
}

/// \brief Standard end() function to make this into a const range
auto ssaps_and_prcs_of_query::end() const -> const_iterator{
	return std::cend( ssap_and_prc_entries );
}

/// \brief TODOCUMENT
optional<reference_wrapper<const ssap_and_prc>> cath::homcheck::best_magic_function(const ssaps_and_prcs_of_query &arg_ssaps_and_prcs ///< TODOCUMENT
                                                                                    ) {
	const auto max_itr = max_element(
		arg_ssaps_and_prcs,
		[&] (const ssap_and_prc &x, const ssap_and_prc &y) {
			return x.get_magic_function_score() < y.get_magic_function_score();
		}
	);
	return ( max_itr != std::cend( arg_ssaps_and_prcs ) ) ? make_optional( cref( *max_itr ) )
	                                                      : none;
}

/// \brief Build a ssaps_and_prcs_of_query from the specified ssap_scores_entries and prc_scores_entries
///
/// \relates ssaps_and_prcs_of_query
ssaps_and_prcs_of_query cath::homcheck::make_ssaps_and_prcs_of_query(const ssap_scores_entry_vec &arg_ssaps, ///< The SSAPs from which the ssaps_and_prcs_of_query should be built
                                                                     const prc_scores_entry_vec  &arg_prcs   ///< The PRCs from which the ssaps_and_prcs_of_query should be built
                                                                     ) {
	// Sanity check the inputs - step 1: check SSAP results have identical name_1s
	const bool ssap_query_ids_identical = all_of(
		arg_ssaps | adjacented,
		[] (const pair<const ssap_scores_entry &, const ssap_scores_entry &> &x) { return x.first.get_name_1() == x.second.get_name_1(); }
	);
	if ( ! ssap_query_ids_identical ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot construct ssaps_and_prcs_of_query from data because the SSAP results have inconsistent query IDs (first is \""
			+ arg_ssaps.front().get_name_1()
			+ "\")"
		));
	}

	// Sanity check the inputs - step 2: check PRC results have identical name_1s
	const bool prc_query_ids_identical = all_of(
		arg_prcs | adjacented,
		[] (const pair<const prc_scores_entry &, const prc_scores_entry &> &x) { return x.first.get_name_1() == x.second.get_name_1(); }
	);
	if ( ! prc_query_ids_identical ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct ssaps_and_prcs_of_query from data because the PRC results have inconsistent query IDs (first is \""
			+ arg_prcs.front().get_name_1()
			+ "\")"
		));
	}

	// Sanity check the inputs - step 3: check the SSAP name_1 matches the PRC name_1
	if ( ! arg_ssaps.empty() && ! arg_prcs.empty() ) {
		const string &ssaps_name_1 = arg_ssaps.front().get_name_1();
		const string &prcs_name_1  =  arg_prcs.front().get_name_1();
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
		irange( 0_z, arg_prcs.size() ),
		[&] (const size_t &x) {
			const auto &the_prc = arg_prcs[ x ];
			return make_pair( make_pair( the_prc.get_name_1(), the_prc.get_name_2() ), x );
		}
	);

	vector<ssap_and_prc> ssap_and_prc_entries;
	for (const ssap_scores_entry &the_ssap : arg_ssaps) {
		const string &query_id = the_ssap.get_name_1();
		const string &match_id = the_ssap.get_name_2();

		const auto prc_index_itr = prc_index_of_ids.find( make_pair( query_id, match_id ) );
		if ( prc_index_itr != std::cend( prc_index_of_ids ) ) {
			ssap_and_prc_entries.emplace_back(
				the_ssap,
				arg_prcs[ prc_index_itr->second ]
			);
		}
	}

	const auto num_ssaps = arg_ssaps.size();
	const auto num_prcs  = arg_prcs.size();
	const auto num_comb  = ssap_and_prc_entries.size();
	const string query_id_str        = ( ! ssap_and_prc_entries.empty() ) ? ( " (query ID : " + ssap_and_prc_entries.front().get_query_id() + ")" ) : "";
	const string unmatched_ssaps_str = ( num_ssaps > num_comb           ) ? ( std::to_string( num_ssaps - num_comb ) + " unmatched SSAP results" ) : "";
	const string unmatched_prcs_str  = ( num_prcs  > num_comb           ) ? ( std::to_string( num_prcs  - num_comb ) + " unmatched PRC results"  ) : "";
	const string conjuction_str      = ( unmatched_ssaps_str.empty() || unmatched_prcs_str.empty() ) ? "" : " and ";
	if ( ! unmatched_ssaps_str.empty() || ! unmatched_prcs_str.empty() ) {
		BOOST_LOG_TRIVIAL( warning ) << "After parsing ssaps_and_prcs_of_query" << query_id_str << ", was left with " << unmatched_ssaps_str << conjuction_str << unmatched_prcs_str;
	}
	return { ssap_and_prc_entries };
}
