/// \file
/// \brief The superfamily_of_domain class definitions

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

#include "superfamily_of_domain.hpp"

#include <filesystem>
#include <fstream>
#include <regex>
#include <unordered_set>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath::common;
using namespace ::cath::homcheck;
using namespace ::cath::homcheck::detail;
using namespace ::std;

using ::boost::algorithm::contains;
using ::boost::algorithm::is_space;
using ::boost::token_compress_on;
using ::std::filesystem::path;

/// \brief The regular expression used to determine whether a string is a valid CATH superfamily ID
const regex is_valid_superfamily_id::SUPERFAMILY_ID_REGEX{ R"(^\d+\.\d+\.\d+\.\d+$)" };

/// \brief The regular expression used to determine whether a string is a valid CATH superfamily ID
const regex is_valid_cath_node_id::NODE_ID_REGEX{ R"(^\d+(\.\d+){0,3}$)" };

/// \brief The string to use in between the fold of a new superfamily and the ID of the domain for which it's being created
const string superfamily_of_domain::NEW_SF_CORE_STRING = ".new_sf_in_fold_of_";

/// \brief Extract the fold from the specified superfamily ID string
///
/// It's assumed that the input is a valid superfamily ID or is constructed by a previous call to fold_of_superfamily_id()
string cath::homcheck::detail::fold_of_superfamily_id(const string &prm_superfamily_id ///< The superfamily ID from which the fold should be extracted
                                                      ) {
	return regex_replace( prm_superfamily_id, regex{ R"(\.[^\.]+$)" }, "" );
}

/// \brief Return whether the specified superfamily is a new superfamily created in this run (rather than a real, original superfamily)
bool superfamily_of_domain::is_created_sf(const string &prm_superfamily ///< The superfamily to query
                                          ) {
	return contains( prm_superfamily, superfamily_of_domain::NEW_SF_CORE_STRING );
}

/// \brief Ctor from a vector<pair<string, string>> where each pair contains domain ID and the corresponding superfamily ID
superfamily_of_domain::superfamily_of_domain(const str_str_pair_vec &prm_sf_of_dom ///< The domain ID -> superfamily ID data from which this superfamily_of_domain should be constructed
                                             ) : sf_of_dom( common::cbegin( prm_sf_of_dom ), common::cend( prm_sf_of_dom ) ) {
	const is_valid_superfamily_id is_valid_sf_pred{};
	for (const auto &x: sf_of_dom) {
		if ( ! is_valid_sf_pred( x.second ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to construct superfamily_of_domain with invalid superfamily ID "+ x.second));
		}
	}
}

/// \brief Getter for the size (ie number of domains for which superfamily information is stored)
size_t superfamily_of_domain::size() const {
	return sf_of_dom.size();
}

/// \brief Return whether this superfamily_of_domain has superfamily information for the specified domain ID
bool superfamily_of_domain::has_superfamily_of_domain(const string &prm_domain_id ///< The domain ID to query
                                                      ) const {
	return sf_of_dom.count( prm_domain_id ) > 0;
}

/// \brief Return whether this superfamily_of_domain has superfamily information for the specified domain ID
///
/// \pre `this->has_superfamily_of_domain( prm_domain_id )`, else and invalid_argument_exception will be thrown
const string & superfamily_of_domain::get_superfamily_of_domain(const string &prm_domain_id ///< The domain ID to query
                                                                ) const {
	const auto find_itr = sf_of_dom.find( prm_domain_id );
	if ( find_itr == common::cend( sf_of_dom ) ) {
		BOOST_THROW_EXCEPTION(
			invalid_argument_exception("Unable to find any entry in superfamily_of_domain for domain ID \""
			+ prm_domain_id
			+ "\""
		));
	}
	return find_itr->second;
}

/// \brief Return whether the specified domain is in a new superfamily created in this run (rather than a real, original superfamily)
///
/// \pre `this->has_superfamily_of_domain( prm_domain_id )`, else and invalid_argument_exception will be thrown
bool superfamily_of_domain::is_in_created_sf(const string &prm_domain_id ///< The domain ID to query
                                             ) const {
	return is_created_sf( get_superfamily_of_domain( prm_domain_id ) );
}

/// \brief Add a new entry representing a new superfamily for the specified domain in the same fold as the second specified domain
///
/// \pre ` ! this->has_superfamily_of_domain( prm_new_domain_id   )`, else and invalid_argument_exception will be thrown
///          this->has_superfamily_of_domain( prm_match_domain_id )`, else and invalid_argument_exception will be thrown
void superfamily_of_domain::add_domain_in_new_sf_in_fold_of_domain(const string &prm_new_domain_id,  ///< The new domain for which an entry should be created
                                                                   const string &prm_match_domain_id ///< The existing domain in whose fold the new domain's new superfamily should be created
                                                                   ) {
	const string &superfamily_id = get_superfamily_of_domain( prm_match_domain_id );
	const string  fold_id        = detail::fold_of_superfamily_id( superfamily_id );

	if ( has_superfamily_of_domain( prm_new_domain_id) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Cannot add entry for domain with ID \""
			+ prm_new_domain_id
			+ "\" into superfamily_of_domain because there it already has an entry"
		));
	}
	sf_of_dom.emplace( prm_new_domain_id, fold_id + NEW_SF_CORE_STRING + prm_match_domain_id );
}

/// \brief Parse the superfamily_of_domain information from the specified istream
///
/// \relates superfamily_of_domain
superfamily_of_domain cath::homcheck::parse_superfamily_of_domain(istream &prm_sf_of_dom_istream ///< The istream from which the superfamily_of_domain information should be parsed
                                                                  ) {
	const is_valid_superfamily_id is_valid_sf_pred{};
	unordered_set<string> prev_seen_domain_ids;
	str_str_pair_vec line_string_pairs;
	string line_string;
	while ( getline( prm_sf_of_dom_istream, line_string ) ) {
		const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
		if ( ! line_parts.empty() ) {
			if ( line_parts.size() != 2 ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception(
					"Cannot parse superfamily_of_domain data because line \""
					+ line_string
					+ "\" does not contain two parts."
				));
			}
			const string &domain_id      = line_parts[ 0 ];
			const string &superfamily_id = line_parts[ 1 ];

			if ( ! is_valid_sf_pred( superfamily_id ) ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception(
					"Second entry in parsed line from superfamily_of_domain data is \""
					+ superfamily_id
					+ "\", which isn't a valid superfamily ID"
				));
			}
			if ( prev_seen_domain_ids.count( domain_id ) > 0 ) {
				BOOST_THROW_EXCEPTION(runtime_error_exception(
					"Duplicated domain_id \""
					+ domain_id
					+ "\" in superfamily_of_domain data"
				));
			}
			prev_seen_domain_ids.insert( domain_id );

			line_string_pairs.emplace_back( domain_id, superfamily_id );
		}
	}

	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return superfamily_of_domain{ line_string_pairs };
}

/// \brief Parse the superfamily_of_domain information from the specified file
///
/// \relates superfamily_of_domain
superfamily_of_domain cath::homcheck::parse_superfamily_of_domain(const path &prm_sf_of_dom_file ///< The file from which the superfamily_of_domain information should be parsed
                                                                  ) {
	ifstream input_ifstream = open_ifstream( prm_sf_of_dom_file );
	const auto sf_of_dom = parse_superfamily_of_domain( input_ifstream );
	input_ifstream.close();
	return sf_of_dom;
}

/// \brief Parse the superfamily_of_domain information from the specified string
///
/// \relates superfamily_of_domain
superfamily_of_domain cath::homcheck::parse_superfamily_of_domain(const string &prm_sf_of_dom_text ///< The string from which the superfamily_of_domain information should be parsed
                                                                  ) {
	istringstream input_ss{  prm_sf_of_dom_text};
	return parse_superfamily_of_domain( input_ss );
}
